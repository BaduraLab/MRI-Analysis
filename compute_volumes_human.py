#!/usr/bin/python3
"""  Compute volumes of human subjects (adjusted with manual cerebellar and lobular annotations)
"""

__author__ = "Enzo Nio"
__version__ = "1.0.0"
__maintainer__ = "Enzo Nio"

import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import glob
import csv
from scipy import stats
from scipy.stats import ttest_ind
from pathlib import Path
import numpy_indexed as npi
from functions import subjectPath2volumeTable
from functions import strPathList2List
from functions import CrawfordHowell
from statsmodels.stats.multitest import multipletests

# Ignore runtime warnings
import warnings
def fxn():
    warnings.warn("runtime", RuntimeWarning)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()



# Define paths
data_path = os.path.join('Data', 'Human', 'Processed')
reference_path = os.path.join('Data', 'Human', 'Reference')
analysis_path = os.path.join('Data', 'Human', 'Analysis')
input_path_list_list = [glob.glob(os.path.join(data_path, '*', '*annotation_orsuit_thrarg_adjusted_lobular.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_subcortical_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_CerebrA_thrarg_adjusted.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_mask_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_AAN_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_allen_thrarg_adjusted.nii.gz'))]
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
                                     'AAN_reoriented.nii.gz'),
                        os.path.join(reference_path,
                                     'allen',
                                     'annotation_full_custom_reoriented.nii.gz')] ###################################### ADD]
structure_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT_plus.csv'),
                       os.path.join(reference_path, 'subcortical', 'subcortical_plus.csv'),
                       os.path.join(reference_path, 'CerebrA', 'CerebrA_plus.csv'),
                       os.path.join(reference_path, 'CerebrA', 'mask_plus.csv'),
                       os.path.join(reference_path, 'AAN', 'AAN_plus.csv'),
                       os.path.join(reference_path, 'allen', 'allen_plus.csv')]
structure_name_list = ['Lobules-SUIT', 'subcortical', 'CerebrA', 'mask', 'AAN', 'allen']
nIterBootstrap = 10000

# Follows
nAnnotation = len(input_path_list_list)
structure_table_list = list()
for iAnnotation in range(nAnnotation):
    structure_table_list.append(pd.read_csv(structure_path_list[iAnnotation]))
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
    print(f'atlas name = {structure_name}')
    # structure_name = os.path.splitext(os.path.basename(structure_path_list[iAnnotation]))[0].split('_')[0]

    # print(annotation_path)

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
    print(f'atlas name = {structure_name}')
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

# For each structure compute WT maskNormalized mean and calculate percemaskNormalizd as a percentage of WT maskNormalized mean
output_table_all['VolumeWTmean'] = np.nan
output_table_all['VolumeWTmeanNormalized'] = np.nan
for name in np.unique(output_table_all['name']):
    WTmean = np.mean(output_table_all.loc[np.logical_and(output_table_all['subject'] != 'patient', output_table_all['name'] == name), 'VolumeMaskNormalized'])
    output_table_all.loc[output_table_all['name'] == name, 'VolumeWTmean'] = WTmean
    output_table_all.loc[output_table_all['name'] == name, 'VolumeWTmeanNormalized'] = (output_table_all['VolumeMaskNormalized'] / WTmean) * 100

output_table_all.to_csv(os.path.join(analysis_path, 'all_volumes_human.csv'))
output_table_cerebellum = output_table_all[output_table_all['atlas'] == 'Lobules-SUIT']
output_table_cerebellum.to_csv(os.path.join(analysis_path, 'cerebellum_volumes_human_human.csv'))

output_table_sn = output_table_all[np.logical_or(output_table_all['name'] == 'substantia nigra pars compacta',
                                                 output_table_all['name'] == 'substantia nigra pars reticula')]
output_table_sn.to_csv(os.path.join(analysis_path, 'sn_volumes_human.csv'))



## Go through each structure in volume table and summarize statistics per structure
output_table_pername_list = list()
uniq_NameAtlas = np.unique(np.vstack([np.array(output_table_all['name'].astype(str)),
                                      np.array(output_table_all['atlas'].astype(str)),
                                      np.array(output_table_all['VolumeInteger'].astype(str))]).astype(str), axis=1)
nNameAtlas = uniq_NameAtlas.shape[1]
for iNameAtlas in range(nNameAtlas):
    nameStruct = uniq_NameAtlas[0, iNameAtlas]
    atlasStruct = uniq_NameAtlas[1, iNameAtlas]
    volintStruct = int(uniq_NameAtlas[2, iNameAtlas])

    logicalStruct = np.logical_and(np.logical_and(output_table_all['name'] == nameStruct,
                                  output_table_all['atlas'] == atlasStruct),
                   output_table_all['VolumeInteger'] == volintStruct)

    output_table_pername = output_table_all[logicalStruct]
    atlas_name = output_table_pername['atlas'].iloc[0]

    controlVolArray = np.array(output_table_pername[output_table_pername['subject'] != 'patient']['Volume'])
    controlMeanVolume = np.mean(controlVolArray)
    controlSDVolume = np.std(controlVolArray)
    patientVolume = float(output_table_pername[output_table_pername['subject'] == 'patient']['Volume'])
    tStat_CH_vol, df_CH_vol, pVal_CH_vol = CrawfordHowell(case=patientVolume, control=controlVolArray)
    tStat_OS_vol, pVal_OS_vol = stats.ttest_1samp(a=controlVolArray, popmean=patientVolume)

    controlVolArray = np.array(output_table_pername[output_table_pername['subject'] != 'patient']['VolumeNormalized'])
    controlMeanVolumeNorm = np.mean(controlVolArray)
    controlSDVolumeNorm = np.std(controlVolArray)
    patientVolumeNorm = float(output_table_pername[output_table_pername['subject'] == 'patient']['VolumeNormalized'])
    tStat_CH_volNorm, df_CH_volNorm, pVal_CH_volNorm = CrawfordHowell(case=patientVolumeNorm, control=controlVolArray)
    tStat_OS_volNorm, pVal_OS_volNorm = stats.ttest_1samp(a=controlVolArray, popmean=patientVolumeNorm)

    controlVolArray = np.array(output_table_pername[output_table_pername['subject'] != 'patient']['VolumeMaskNormalized'])
    controlMeanVolumeMaskNorm = np.mean(controlVolArray)
    controlSDVolumeMaskNorm = np.std(controlVolArray)
    patientVolumeMaskNorm = float(output_table_pername[output_table_pername['subject'] == 'patient']['VolumeMaskNormalized'])
    tStat_CH_volNormMask, df_CH_volNormMask, pVal_CH_volNormMask = CrawfordHowell(case=patientVolumeMaskNorm, control=controlVolArray)
    tStat_OS_volNormMask, pVal_OS_volNormMask = stats.ttest_1samp(a=controlVolArray, popmean=patientVolumeMaskNorm)

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


    # TODO: implement CrawfordHowell as well as single case
    # CrawfordHowell < - function(case, control)
    # {
    #     tval < - (case - mean(control)) / (sd(control) * sqrt((length(control) + 1) / length(control)))
    # degfree < - length(control) - 1
    # pval < - 2 * (1 - pt(abs(tval), df=degfree))  # two-tailed p-value
    # result < - data.frame(t=tval, df=degfree, p=pval)
    # return (result)
    # }
    # tStat_CH, df_CH, pVal_CH = CrawfordHowell(case=patientVolume, control=controlVolArray)
    # tStat_OS, pVal_OS = stats.ttest_1samp(a=controlVolArray, popmean=patientVolume)





    VolumeInteger_val = np.array(output_table_pername['VolumeInteger'])[0]

    output_table_pername_list.append(pd.DataFrame({'name': [nameStruct],
                                                   'atlas': [atlasStruct],
                                                   'VolumeInteger': [VolumeInteger_val],
                                                   'pVal': [pVal_CH_vol],
                                                   'cohenD': [cohenD],
                                                   'cohenDMaskNorm': [cohenDMaskNorm],
                                                   'pValMaskNorm': [pVal_CH_volNormMask],
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
                                                   'patientVolumeMaskNorm': [patientVolumeMaskNorm],
                                                   'tStat_CH_vol': [tStat_CH_vol],
                                                   'pVal_CH_vol': [pVal_CH_vol],
                                                   'tStat_OS_vol': [tStat_OS_vol],
                                                   'pVal_OS_vol': [pVal_OS_vol],
                                                   'tStat_CH_volNorm': [tStat_CH_volNorm],
                                                   'pVal_CH_volNorm': [pVal_CH_volNorm],
                                                   'tStat_OS_volNorm': [tStat_OS_volNorm],
                                                   'pVal_OS_volNorm': [pVal_OS_volNorm],
                                                   'tStat_CH_volNormMask': [tStat_CH_volNormMask],
                                                   'pVal_CH_volNormMask': [pVal_CH_volNormMask],
                                                   'tStat_OS_volNormMask': [tStat_OS_volNormMask],
                                                   'pVal_OS_volNormMask': [pVal_OS_volNormMask]}))


output_table_pername = pd.concat(output_table_pername_list, ignore_index=True)
output_table_pername = output_table_pername.sort_values(by='cohenDMaskNorm', ascending=False)
output_table_pername['relFracDiff'] = output_table_pername['fractionDifference']/float(output_table_pername[output_table_pername['name']=='mask']['fractionDifference'])
output_table_pername.loc[np.logical_not(np.isnan(output_table_pername['pVal'])), 'pValFDR'] = \
    multipletests(output_table_pername.loc[np.logical_not(np.isnan(output_table_pername['pVal'])), 'pVal'], method='fdr_bh')[1]
output_table_pername.loc[np.logical_not(np.isnan(output_table_pername['pValMaskNorm'])), 'pValMaskNormFDR'] = \
    multipletests(output_table_pername.loc[np.logical_not(np.isnan(output_table_pername['pValMaskNorm'])), 'pValMaskNorm'], method='fdr_bh')[1]
# output_table_pername = output_table_pername[['name', 'atlas', 'VolumeInteger', 'pVal', 'pValFDR', 'cohenD',
#                                              'cohenDMaskNorm', 'pValMaskNorm', 'cohenD_CI', 'cohenDMaskNorm_CI', 'fractionDifference', 'controlMeanVolume']]
output_table_pername.insert(4, 'pValFDR_temp', output_table_pername['pValFDR'])
output_table_pername.drop('pValFDR', axis=1, inplace=True)
output_table_pername.rename(columns={'pValFDR_temp': 'pValFDR'}, inplace=True)
output_table_pername.insert(8, 'pValMaskNormFDR_temp', output_table_pername['pValMaskNormFDR'])
output_table_pername.drop('pValMaskNormFDR', axis=1, inplace=True)
output_table_pername.rename(columns={'pValMaskNormFDR_temp': 'pValMaskNormFDR'}, inplace=True)
output_table_pername.to_csv(os.path.join(analysis_path, 'pername_volumes_human.csv'))

output_table_pername_cerebellum = output_table_pername[output_table_pername['atlas'] == 'Lobules-SUIT']
output_table_pername_cerebellum.loc[np.logical_not(np.isnan(output_table_pername_cerebellum['pVal'])), 'pValFDR'] = \
    multipletests(output_table_pername_cerebellum.loc[np.logical_not(np.isnan(output_table_pername_cerebellum['pVal'])), 'pVal'], method='fdr_bh')[1]
output_table_pername_cerebellum.loc[np.logical_not(np.isnan(output_table_pername_cerebellum['pValMaskNorm'])), 'pValMaskNormFDR'] = \
    multipletests(output_table_pername_cerebellum.loc[np.logical_not(np.isnan(output_table_pername_cerebellum['pValMaskNorm'])), 'pValMaskNorm'], method='fdr_bh')[1]
output_table_pername_cerebellum.to_csv(os.path.join(analysis_path, 'pername_cerebellum_volumes_human.csv'))

output_table_pername_subcortical = output_table_pername[output_table_pername['atlas'] == 'subcortical']
output_table_pername_subcortical.loc[np.logical_not(np.isnan(output_table_pername_subcortical['pVal'])), 'pValFDR'] = \
    multipletests(output_table_pername_subcortical.loc[np.logical_not(np.isnan(output_table_pername_subcortical['pVal'])), 'pVal'], method='fdr_bh')[1]
output_table_pername_subcortical.loc[np.logical_not(np.isnan(output_table_pername_subcortical['pValMaskNorm'])), 'pValMaskNormFDR'] = \
    multipletests(output_table_pername_subcortical.loc[np.logical_not(np.isnan(output_table_pername_subcortical['pValMaskNorm'])), 'pValMaskNorm'], method='fdr_bh')[1]
output_table_pername_subcortical.to_csv(os.path.join(analysis_path, 'pername_subcortical_volumes_human.csv'))



# Create 3 images based on pername for each atlas: 1) Map by abs CohenD 2) Same but increased (CohenD<0) 3) Same but decreased (CohenD>0) 
# Create 3 additional images but for CohenDMaskNorm
for iAnnotation in range(nAnnotation):
    annotation_path = annotation_path_list[iAnnotation]
    annotation_image = nib.load(annotation_path)
    annotation = annotation_image.get_fdata()
    structure_table = structure_table_list[iAnnotation]
    structure_name = structure_name_list[iAnnotation]
    print(f'atlas name = {structure_name}')
    output_table_pername_peratlas = output_table_pername[output_table_pername['atlas'] == structure_name]

    map_from = np.array(output_table_pername_peratlas['VolumeInteger'])
    
    for j in range(2):
        if j==0:
            mapArr = np.array(output_table_pername_peratlas['cohenD'])
        else:
            mapArr = np.array(output_table_pername_peratlas['cohenDMaskNorm'])
            
        for i in range(3):
            print(i)
    
            annotation_CD_path = os.path.join(analysis_path, annotation_path.split(os.sep)[-1].split('.')[0])
            if i == 0:
                map_to = np.abs(mapArr)
                annotation_CD_path = annotation_CD_path + '_CD' + '.nii.gz'
            elif i == 1:
                map_to = np.abs(mapArr)*(mapArr < 0)
                annotation_CD_path = annotation_CD_path + '_CD_volIncrease' + '.nii.gz'
            elif i == 2:
                map_to = np.abs(mapArr)*(mapArr > 0)
                annotation_CD_path = annotation_CD_path + '_CD_volDecrease' + '.nii.gz'
    
            annotation_remapped = np.round(annotation) # always annotation so should never be non-integer
            # input = input.astype(int) # always annotation so should never be non-integer
            annotation_remapped_shape = annotation_remapped.shape
            annotation_remapped = annotation_remapped.reshape(-1)
            annotation_remapped = npi.remap(annotation_remapped, map_from, map_to)
            annotation_remapped = annotation_remapped.reshape(annotation_remapped_shape)
            # annotation_remapped = remap_3D(annotation, map_from_filtered.astype(int), map_to_filtered)
    
            output_image = nib.Nifti1Image(annotation_remapped,
                                           annotation_image.affine)
            nib.save(output_image, annotation_CD_path)



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

# #### #### #### ########################################################################
# ## Go through each structure in volume table and summarize statistics per structure
# uniq_NameAtlas_pername = np.unique(np.vstack([np.array(output_table_pername['name'].astype(str)),
#                                               np.array(output_table_pername['atlas'].astype(str)),
#                                               np.array(output_table_pername['VolumeInteger'].astype(str))]).astype(str), axis=1)
# nNameAtlas_pername = uniq_NameAtlas_pername.shape[1]
# for iNameAtlas in range(nNameAtlas_pername):
#     nameStruct = uniq_NameAtlas_pername[0, iNameAtlas]
#     atlasStruct = uniq_NameAtlas_pername[1, iNameAtlas]
#     volintStruct = int(uniq_NameAtlas[2, iNameAtlas])
# 
#     logicalStruct = np.logical_and(np.logical_and(output_table_all['name'] == nameStruct,
#                                   output_table_all['atlas'] == atlasStruct),
#                    output_table_all['VolumeInteger'] == volintStruct)
# 
#     logicalStruct_pername = np.logical_and(np.logical_and(output_table_pername['name'] == nameStruct,
#                                   output_table_pername['atlas'] == atlasStruct),
#                    output_table_pername['VolumeInteger'] == volintStruct)
# 
#     cohenD_current = float(output_table_pername.loc[logicalStruct_pername, 'cohenD'])
#     cohenDMaskNorm_current = float(output_table_pername.loc[logicalStruct_pername, 'cohenDMaskNorm'])
# 
#     fig = plt.figure()
#     output_table_all_plot = output_table_all[logicalStruct]
#     ax = output_table_all_plot.plot.bar(x='subject',
#                                         y='VolumeNormalized',
#                                         rot=0,
#                                         color=[[0.4, 0.4, 0.85],
#                                                [0.4, 0.4, 0.85],
#                                                [0.4, 0.4, 0.85],
#                                                [0.4, 0.4, 0.85],
#                                                [0.4, 0.4, 0.85],
#                                                [0.85, 0.4, 0.4]])
#     plt.ylabel('Volume Normalized to Reference')
#     plt.title(nameStruct + '_' + atlasStruct + ', CohenD = ' + format('%.2f'%cohenD_current))
#     ax.get_legend().remove()
#     # plt.show()
#     plt.savefig(os.path.join(analysis_path,
#                              'barplots',
#                              'CohenD_'+ format('%.2f'%cohenD_current) + '_' + nameStruct + '_' + atlasStruct + '_volume_barplot.png'))
#     plt.close('all')
# 
# 
#     fig = plt.figure()
#     output_table_all_plot = output_table_all[logicalStruct]
#     ax = output_table_all_plot.plot.bar(x='subject',
#                                         y='VolumeMaskNormalized',
#                                         rot=0,
#                                         color=[[0.4, 0.4, 0.85],
#                                                [0.4, 0.4, 0.85],
#                                                [0.4, 0.4, 0.85],
#                                                [0.4, 0.4, 0.85],
#                                                [0.4, 0.4, 0.85],
#                                                [0.85, 0.4, 0.4]])
#     plt.ylabel('Volume Percentage')
#     plt.title(nameStruct + '_' + atlasStruct + ', CohenD = ' + format('%.2f'%cohenDMaskNorm_current))
#     ax.get_legend().remove()
#     # plt.show()
#     plt.savefig(os.path.join(analysis_path,
#                              'barplots',
#                              'CohenD_'+ format('%.2f'%cohenDMaskNorm_current) + '_' + nameStruct + '_' + atlasStruct + '_volume_barplot_MaskNorm.png'))
#     plt.close('all')
# 
# 
#     #
# 
#     fig = plt.figure()
#     output_table_all_plot = output_table_all[logicalStruct]
#     output_table_all_plot.loc[np.array(output_table_all_plot['subject'] == 'patient'), 'Genotype'] = 'patient'
#     output_table_all_plot.loc[np.array(output_table_all_plot['subject'] != 'patient'), 'Genotype'] = 'control'
#     ax = output_table_all_plot[['VolumeNormalized', 'Genotype']].boxplot(by=['Genotype'])
#     plt.ylabel('Volume Normalized to Reference')
#     plt.title(nameStruct + '_' + atlasStruct + ', CohenD = ' + format('%.2f'%cohenD_current))
#     plt.suptitle('')
#     # ax.get_legend().remove()
#     # plt.show()
#     # ax.set_aspect(1.5)
#     plt.savefig(os.path.join(analysis_path,
#                              'boxplots',
#                              'CohenD_' + format('%.2f'%cohenD_current) + '_' + nameStruct + '_' + atlasStruct + '_volume_boxplot.png'))
#     plt.close('all')
# 
# 
# 
#     fig = plt.figure()
#     output_table_all_plot = output_table_all[logicalStruct]
#     output_table_all_plot.loc[np.array(output_table_all_plot['subject'] == 'patient'), 'Genotype'] = 'patient'
#     output_table_all_plot.loc[np.array(output_table_all_plot['subject'] != 'patient'), 'Genotype'] = 'control'
#     ax = output_table_all_plot[['VolumeMaskNormalized', 'Genotype']].boxplot(by=['Genotype'])
#     plt.ylabel('Volume Percentage')
#     plt.title(nameStruct + '_' + atlasStruct + ', CohenD = ' + format('%.2f'%cohenDMaskNorm_current))
#     plt.suptitle('')
#     # ax.get_legend().remove()
#     # plt.show()
#     # ax.set_aspect(1.5)
#     plt.savefig(os.path.join(analysis_path,
#                              'boxplots',
#                              'CohenD_'+ format('%.2f'%cohenDMaskNorm_current) + '_' + nameStruct + '_' + atlasStruct + '_volume_boxplot_MaskNorm.png'))
#     plt.close('all')
