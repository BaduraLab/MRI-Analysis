import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import glob
import csv
from scipy.stats import ttest_ind
from pathlib import Path
import numpy_indexed as npi
from functions import subjectPath2volumeTable
from functions import strPathList2List



# Define paths
data_path = os.path.join('Data', 'Human', 'Processed')
reference_path = os.path.join('Data', 'Human', 'Reference')
analysis_path = os.path.join('Data', 'Human', 'Analysis')
input_path_list_list = [glob.glob(os.path.join(data_path, '*', '*annotation_orsuit_thrarg_adjusted_lobular.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_subcortical_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_CerebrA_thrarg_adjusted.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_mask_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_AAN_thrarg.nii.gz'))]
annotation_path_list = [os.path.join(reference_path,
                                     'suit',
                                     'atlasesSUIT',
                                     'Lobules-SUIT.nii'),
                        os.path.join(reference_path,
                                     'subcortical',
                                     'prob_atlas_bilateral_thrarg_0.4.nii.gz'),
                        os.path.join(reference_path,
                                     'CerebrA',
                                     'mni_icbm152_CerebrA_tal_nlin_sym_09c_reoriented.nii.gz'),
                        os.path.join(reference_path,
                                     'CerebrA',
                                     'mni_icbm152_t1_tal_nlin_sym_09c_mask_reoriented.nii.gz'),
                        os.path.join(reference_path,
                                     'AAN',
                                     'AAN_reoriented.nii.gz')] ###################################### ADD]
structure_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT_plus.csv'),
                       os.path.join(reference_path, 'subcortical', 'subcortical_plus.csv'),
                       os.path.join(reference_path, 'CerebrA', 'CerebrA_plus.csv'),
                       os.path.join(reference_path, 'CerebrA', 'mask_plus.csv'),
                       os.path.join(reference_path, 'AAN', 'AAN_plus.csv')]
structure_name_list = ['Lobules-SUIT', 'subcortical', 'CerebrA', 'mask', 'AAN']
nIterBootstrap = 10000

# Follows
nAnnotation = len(input_path_list_list)
structure_table_list = [pd.read_csv(structure_path_list[0]),
                        pd.read_csv(structure_path_list[1]),
                        pd.read_csv(structure_path_list[2]),
                        pd.read_csv(structure_path_list[3])]
pd.concat(structure_table_list, sort=False).to_csv(os.path.join(reference_path, 'structure_table.csv'))

# for iAnnotation in range(3):
#     ind = iAnnotation + 1
#     structure_table_list[ind]['VolumeIntegerPath'] = \
#         [[iVol] for iVol in structure_table_list[ind]['VolumeInteger']]
#     structure_table_list[ind].to_csv(structure_path_list[ind].split('.')[0]+'_plus.csv')


## Calculate volumes for all atlas files
output_table_list = list()
for iAnnotation in range(nAnnotation):
    annotation_path = annotation_path_list[iAnnotation]
    structure_table = structure_table_list[iAnnotation]
    structure_name = structure_name_list[iAnnotation]
    # structure_name = os.path.splitext(os.path.basename(structure_path_list[iAnnotation]))[0].split('_')[0]

    print(annotation_path)

    # Compute volumes
    output_table = subjectPath2volumeTable(annotation_path)

    # Add structure names to output table by using common volume integers
    output_table = pd.merge(left=output_table,
                            right=structure_table_list[iAnnotation].loc[:, ['VolumeInteger', 'name', 'VolumeIntegerPath']],
                            left_on='VolumeInteger',
                            right_on='VolumeInteger',
                            how='right')

    # Any volumes with NaN values should be filled in with volumes that contain that volume as a parent
    non_nanVol = np.logical_not(np.isnan(output_table['Volume']))
    for iVol in range(output_table.shape[0]):
        if np.isnan(output_table.iloc[iVol]['Volume']):
            isVol1_Child = list()
            for iVol2 in range(output_table.shape[0]):
                strPathList = strPathList2List(output_table.iloc[iVol2]['VolumeIntegerPath'])[:-1]
                VolumeInteger_iVol = output_table.iloc[iVol]['VolumeInteger']
                isVol1_Child.append(np.any(np.isin(strPathList, VolumeInteger_iVol)))
            output_table.loc[iVol, 'VoxelNumber'] = np.sum(output_table.loc[np.logical_and(isVol1_Child, non_nanVol), 'VoxelNumber'])
            output_table.loc[iVol, 'Volume'] = np.sum(output_table.loc[np.logical_and(isVol1_Child, non_nanVol), 'Volume'])

    # Add metadata
    output_table['atlas'] = structure_name

    # Add table to table list
    output_table_list.append(output_table)

output_table_all = pd.concat(output_table_list, ignore_index=True)
output_table_all.to_csv(os.path.join(analysis_path, 'reference_volumes_human.csv'))
reference_table = output_table_all.copy()
reference_cerebellum_table = reference_table[reference_table['atlas'] == 'Lobules-SUIT']
reference_cerebellum_table.to_csv(os.path.join(analysis_path, 'reference_cerebellum_volumes_human.csv'))
# reference_table[reference_table['']=='']



## Calculate volumes for all input files
output_table_list = list()
for iAnnotation in range(nAnnotation):
    input_path_list = input_path_list_list[iAnnotation]
    structure_table = structure_table_list[iAnnotation]
    structure_name = structure_name_list[iAnnotation]
    # structure_name = os.path.splitext(os.path.basename(structure_path_list[iAnnotation]))[0]

    for iInput, Input in enumerate(input_path_list):
        input_name = Input.split(os.sep)[-2]
        print('subject '+str(iInput)+': '+input_name)

        # Compute volumes
        output_table = subjectPath2volumeTable(Input)

        # Add structure names to output table by using common volume integers
        output_table = pd.merge(left=output_table,
                                right=structure_table_list[iAnnotation].loc[:,
                                      ['VolumeInteger', 'name', 'VolumeIntegerPath']],
                                left_on='VolumeInteger',
                                right_on='VolumeInteger',
                                how='right')

        # Any volumes with NaN values should be filled in with volumes that contain that volume as a parent
        non_nanVol = np.logical_not(np.isnan(output_table['Volume']))
        for iVol in range(output_table.shape[0]):
            if np.isnan(output_table.iloc[iVol]['Volume']):
                isVol1_Child = list()
                for iVol2 in range(output_table.shape[0]):
                    strPathList = strPathList2List(output_table.iloc[iVol2]['VolumeIntegerPath'])[:-1]
                    VolumeInteger_iVol = output_table.iloc[iVol]['VolumeInteger']
                    isVol1_Child.append(np.any(np.isin(strPathList, VolumeInteger_iVol)))
                output_table.loc[iVol, 'VoxelNumber'] = np.sum(
                    output_table.loc[np.logical_and(isVol1_Child, non_nanVol), 'VoxelNumber'])
                output_table.loc[iVol, 'Volume'] = np.sum(
                    output_table.loc[np.logical_and(isVol1_Child, non_nanVol), 'Volume'])

        # Add metadata
        output_table['subject'] = input_name
        output_table['atlas'] = structure_name

        # Add table to table list
        output_table_list.append(output_table)

output_table_all = pd.concat(output_table_list, ignore_index=True)

reference_table['rVolume'] = reference_table['Volume']
output_table_all = pd.merge(left=output_table_all,
                            right=reference_table.loc[:, ['atlas', 'VolumeInteger', 'name', 'rVolume']],
                            left_on=['atlas', 'VolumeInteger', 'name'],
                            right_on=['atlas', 'VolumeInteger', 'name'])
output_table_all['VolumeNormalized'] = output_table_all['Volume']/output_table_all['rVolume']

# For each subject add mask volumes for  normalization
output_table_all['maskVolume'] = np.nan
for subject in np.unique(output_table_all['subject']):
    output_table_all.loc[output_table_all['subject'] == subject, 'maskVolume'] = float(output_table_all.loc[np.logical_and(output_table_all['subject'] == subject, output_table_all['atlas'] == 'mask'), 'Volume'])

output_table_all['VolumeMaskNormalized'] = output_table_all['Volume']/output_table_all['maskVolume']
output_table_all.to_csv(os.path.join(analysis_path, 'all_volumes_human.csv'))
output_table_cerebellum = output_table_all[output_table_all['atlas'] == 'Lobules-SUIT']
output_table_cerebellum.to_csv(os.path.join(analysis_path, 'cerebellum_volumes_human_human.csv'))

output_table_sn = output_table_all[np.logical_or(output_table_all['name'] == 'substantia nigra pars compacta',
                                                 output_table_all['name'] == 'substantia nigra pars reticula')]
output_table_sn.to_csv(os.path.join(analysis_path, 'sn_volumes_human.csv'))



## Go through each structure in volume table and summarize statistics per structure
output_table_pername_list = list()
for nameStruct in np.unique(np.array(output_table_all['name'].astype('category'))):

    output_table_pername = output_table_all[output_table_all['name'] == nameStruct]
    atlas_name = output_table_pername['atlas'].iloc[0]

    controlVolArray = np.array(output_table_pername[output_table_pername['subject'] != 'patient']['Volume'])
    controlMeanVolume = np.mean(controlVolArray)
    controlSDVolume = np.std(controlVolArray)
    patientVolume = float(output_table_pername[output_table_pername['subject'] == 'patient']['Volume'])

    controlVolArray = np.array(output_table_pername[output_table_pername['subject'] != 'patient']['VolumeNormalized'])
    controlMeanVolumeNorm = np.mean(controlVolArray)
    controlSDVolumeNorm = np.std(controlVolArray)
    patientVolumeNorm = float(output_table_pername[output_table_pername['subject'] == 'patient']['VolumeNormalized'])

    controlVolArray = np.array(output_table_pername[output_table_pername['subject'] != 'patient']['VolumeMaskNormalized'])
    controlMeanVolumeMaskNorm = np.mean(controlVolArray)
    controlSDVolumeMaskNorm = np.std(controlVolArray)
    patientVolumeMaskNorm = float(output_table_pername[output_table_pername['subject'] == 'patient']['VolumeMaskNormalized'])

    fractionDifferenceVolume = (patientVolume - controlMeanVolume) / controlMeanVolume
    nControl = len(controlVolArray)
    N = nControl + 1
    hedge_correction = (N - 3) / (N - 2.25)
    cohenD = ((controlMeanVolume - patientVolume) / controlSDVolume) * hedge_correction
    # print(f'nControl = {nControl}, N = {N}, hedge_correction = {hedge_correction}, cohenD = {cohenD}')

    cohenD_BS = np.empty(nIterBootstrap)
    for iBS in range(nIterBootstrap):
        controlVolArrayNormal = np.random.normal(controlMeanVolume, controlSDVolume, nControl)
        controlMeanVolumeBS = np.mean(controlVolArrayNormal)
        controlSDVolumeBS = np.std(controlVolArrayNormal)
        cohenD_BS[iBS] = ((controlMeanVolumeBS - patientVolume) / controlSDVolumeBS) * hedge_correction
    cohenD_CI = [np.quantile(cohenD_BS, .025), np.quantile(cohenD_BS, .975)]
    # cohenD_check = [np.quantile(cohenD_BS, .5), np.median(cohenD_BS)]
    # print(cohenD_check)

    cohenDMaskNorm = ((controlMeanVolumeMaskNorm - patientVolumeMaskNorm) / controlSDVolumeMaskNorm) * hedge_correction

    cohenD_BS = np.empty(nIterBootstrap)
    for iBS in range(nIterBootstrap):
        controlVolArrayNormalMaskNorm = np.random.normal(controlMeanVolumeMaskNorm, controlSDVolumeMaskNorm, nControl)
        controlMeanVolumeBS = np.mean(controlVolArrayNormalMaskNorm)
        controlSDVolumeBS = np.std(controlVolArrayNormalMaskNorm)
        cohenD_BS[iBS] = ((controlMeanVolumeBS - patientVolumeMaskNorm) / controlSDVolumeBS) * hedge_correction
    cohenDMaskNorm_CI = [np.quantile(cohenD_BS, .025), np.quantile(cohenD_BS, .975)]



    VolumeInteger_val = np.array(output_table_pername['VolumeInteger'])[0]

    output_table_pername_list.append(pd.DataFrame({'name': [nameStruct],
                                                  'atlas': [atlas_name],
                                                   'VolumeInteger': [VolumeInteger_val],
                                                  'cohenD': [cohenD],
                                                  'cohenDMaskNorm': [cohenDMaskNorm],
                                                  'cohenD_CI': [cohenD_CI],
                                                  'cohenDMaskNorm_CI': [cohenDMaskNorm_CI],
                                                  'fractionDifference': [fractionDifferenceVolume],
                                                  'controlMeanVolume': [controlMeanVolume],
                                                  'controlSDVolume': [controlSDVolume],
                                                  'patientVolume': [patientVolume],
                                                  'controlMeanVolumeNorm': [controlMeanVolumeNorm],
                                                  'controlSDVolumeNorm': [controlSDVolumeNorm],
                                                  'patientVolumeNorm': [patientVolumeNorm],
                                                  'controlMeanVolumeMaskNorm': [controlMeanVolumeMaskNorm],
                                                  'controlSDVolumeMaskNorm': [controlSDVolumeMaskNorm],
                                                  'patientVolumeMaskNorm': [patientVolumeMaskNorm]}))

output_table_pername = pd.concat(output_table_pername_list, ignore_index=True)
output_table_pername = output_table_pername.sort_values(by='cohenDMaskNorm', ascending=False)
output_table_pername['relFracDiff'] = output_table_pername['fractionDifference']/float(output_table_pername[output_table_pername['name']=='mask']['fractionDifference'])
output_table_pername.to_csv(os.path.join(analysis_path, 'pername'+'_volumes_human.csv'))
output_table_pername_cerebellum = output_table_pername[output_table_pername['atlas'] == 'Lobules-SUIT']
output_table_pername_cerebellum.to_csv(os.path.join(analysis_path, 'pername_cerebellum_volumes_human.csv'))



# Create 3 images based on pername for each atlas: 1) Map by abs relFracDiff 2) Same but increased (relFracDiff>0) 3) Same but decreased (relFracDiff<0)
for iAnnotation in range(nAnnotation):
    print(iAnnotation)
    annotation_path = annotation_path_list[iAnnotation]
    annotation_image = nib.load(annotation_path)
    annotation = annotation_image.get_fdata()
    structure_table = structure_table_list[iAnnotation]
    structure_name = structure_name_list[iAnnotation]
    output_table_pername_peratlas = output_table_pername[output_table_pername['atlas'] == structure_name]

    map_from = np.array(output_table_pername_peratlas['VolumeInteger'])

    for i in range(3):
        print(i)

        annotation_aFD_path = os.path.join(analysis_path, annotation_path.split(os.sep)[-1].split('.')[0])
        if i == 0:
            map_to = np.abs(np.array(output_table_pername_peratlas['fractionDifference']))
            annotation_aFD_path = annotation_aFD_path + '_aFD' + '.nii.gz'
        elif i == 1:
            map_to = np.abs(np.array(output_table_pername_peratlas['fractionDifference']))*(np.array(output_table_pername_peratlas['fractionDifference']) > 0)
            annotation_aFD_path = annotation_aFD_path + '_aFD_volIncrease' + '.nii.gz'
        elif i == 2:
            map_to = np.abs(np.array(output_table_pername_peratlas['fractionDifference']))*(np.array(output_table_pername_peratlas['fractionDifference']) < 0)
            annotation_aFD_path = annotation_aFD_path + '_aFD_volDecrease' + '.nii.gz'

        annotation_remapped = np.round(annotation) # always annotation so should never be non-integer
        # input = input.astype(int) # always annotation so should never be non-integer
        annotation_remapped_shape = annotation_remapped.shape
        annotation_remapped = annotation_remapped.reshape(-1)
        annotation_remapped = npi.remap(annotation_remapped, map_from, map_to)
        annotation_remapped = annotation_remapped.reshape(annotation_remapped_shape)
        # annotation_remapped = remap_3D(annotation, map_from_filtered.astype(int), map_to_filtered)

        output_image = nib.Nifti1Image(annotation_remapped,
                                       annotation_image.affine)
        nib.save(output_image, annotation_aFD_path)



## Plotting
# cmap = matplotlib.cm.get_cmap('lines')
# VOIs = ['substantia nigra pars reticula',
#         'substantia nigra pars compacta',
#         'nucleus accumbens',
#         'mask',
#         'Right_VIIb',
#         'Left_VIIb',
#         'Right_VIIIa',
#         'Right_VIIIb',
#         'Left_VIIIa',
#         'Left_VIIIb']
# VOIs = ['substantia nigra pars reticula',
#         'substantia nigra pars compacta',
#         'nucleus accumbens',
#         'mask',
#         'Left_CrusI',
#         'Right_CrusI',
#         'Vermis_CrusI',
#         'Left_CrusII',
#         'Right_CrusII',
#         'Left_X',
#         'Right_X',
#         'Left_VIIb',
#         'Right_VIIIa']
VOIs = list(output_table_pername['name'])
for iVOI in range(len(VOIs)):
    cohenD_current = float(output_table_pername.loc[output_table_pername['name'] == VOIs[iVOI], 'cohenD'])
    cohenDMaskNorm_current = float(output_table_pername.loc[output_table_pername['name'] == VOIs[iVOI], 'cohenDMaskNorm'])

    output_table_all_plot = output_table_all[output_table_all['name'] == VOIs[iVOI]]
    ax = output_table_all_plot.plot.bar(x='subject',
                                        y='VolumeNormalized',
                                        rot=0,
                                        color=[[0.4, 0.4, 0.85],
                                               [0.4, 0.4, 0.85],
                                               [0.4, 0.4, 0.85],
                                               [0.4, 0.4, 0.85],
                                               [0.4, 0.4, 0.85],
                                               [0.85, 0.4, 0.4]])
    plt.ylabel('Volume Normalized to Reference')
    plt.title(VOIs[iVOI]+', CohenD = ' + format('%.2f'%cohenD_current))
    ax.get_legend().remove()
    # plt.show()
    plt.savefig(os.path.join(analysis_path,
                             'barplots',
                             'CohenD_'+ format('%.2f'%cohenD_current) +'_'+VOIs[iVOI]+'_volume_barplot.png'))


    output_table_all_plot = output_table_all[output_table_all['name'] == VOIs[iVOI]]
    ax = output_table_all_plot.plot.bar(x='subject',
                                        y='VolumeMaskNormalized',
                                        rot=0,
                                        color=[[0.4, 0.4, 0.85],
                                               [0.4, 0.4, 0.85],
                                               [0.4, 0.4, 0.85],
                                               [0.4, 0.4, 0.85],
                                               [0.4, 0.4, 0.85],
                                               [0.85, 0.4, 0.4]])
    plt.ylabel('Volume Percentage')
    plt.title(VOIs[iVOI]+', CohenD = ' + format('%.2f'%cohenDMaskNorm_current))
    ax.get_legend().remove()
    # plt.show()
    plt.savefig(os.path.join(analysis_path,
                             'barplots',
                             'CohenD_'+ format('%.2f'%cohenDMaskNorm_current) +'_'+VOIs[iVOI]+'_volume_barplot_MaskNorm.png'))



    #


    output_table_all_plot = output_table_all[output_table_all['name'] == VOIs[iVOI]]
    output_table_all_plot.loc[np.array(output_table_all_plot['subject'] == 'patient'), 'Genotype'] = 'patient'
    output_table_all_plot.loc[np.array(output_table_all_plot['subject'] != 'patient'), 'Genotype'] = 'control'
    ax = output_table_all_plot[['VolumeNormalized', 'Genotype']].boxplot(by=['Genotype'])
    plt.ylabel('Volume Normalized to Reference')
    plt.title(VOIs[iVOI]+', CohenD = ' + format('%.2f'%cohenD_current))
    plt.suptitle('')
    # ax.get_legend().remove()
    # plt.show()
    # ax.set_aspect(1.5)
    plt.savefig(os.path.join(analysis_path,
                             'boxplots',
                             'CohenD_' + format('%.2f'%cohenD_current) + '_' + VOIs[iVOI] + '_volume_boxplot.png'))


    output_table_all_plot = output_table_all[output_table_all['name'] == VOIs[iVOI]]
    output_table_all_plot.loc[np.array(output_table_all_plot['subject'] == 'patient'), 'Genotype'] = 'patient'
    output_table_all_plot.loc[np.array(output_table_all_plot['subject'] != 'patient'), 'Genotype'] = 'control'
    ax = output_table_all_plot[['VolumeMaskNormalized', 'Genotype']].boxplot(by=['Genotype'])
    plt.ylabel('Volume Percentage')
    plt.title(VOIs[iVOI]+', CohenD = ' + format('%.2f'%cohenDMaskNorm_current))
    plt.suptitle('')
    # ax.get_legend().remove()
    # plt.show()
    # ax.set_aspect(1.5)
    plt.savefig(os.path.join(analysis_path,
                             'boxplots',
                             'CohenD_'+ format('%.2f'%cohenDMaskNorm_current) +'_'+VOIs[iVOI]+'_volume_boxplot_MaskNorm.png'))
