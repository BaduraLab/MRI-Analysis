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



# Define paths
data_path = os.path.join('Data', 'Human', 'Processed')
reference_path = os.path.join('Data', 'Human', 'Reference')
analysis_path = os.path.join('Data', 'Human', 'Analysis')
input_path_list_list = [glob.glob(os.path.join(data_path, '*', '*annotation_orsuit_thrarg_adjusted_lobular.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_subcortical_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_CerebrA_thrarg_adjusted.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_mask_thrarg.nii.gz'))]
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
                                     'mni_icbm152_t1_tal_nlin_sym_09c_mask_reoriented.nii.gz')]
structure_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT.csv'),
                       os.path.join(reference_path, 'subcortical', 'subcortical.csv'),
                       os.path.join(reference_path, 'CerebrA', 'CerebrA.csv'),
                       os.path.join(reference_path, 'CerebrA', 'mask.csv')]
structure_name_list = ['Lobules-SUIT', 'subcortical', 'CerebrA', 'mask']

# Follows
nAnnotation = len(input_path_list_list)
structure_table_list = [pd.read_csv(structure_path_list[0]),
                        pd.read_csv(structure_path_list[1]),
                        pd.read_csv(structure_path_list[2]),
                        pd.read_csv(structure_path_list[3])]



## Calculate volumes for all atlas files
output_table_list = list()
for iAnnotation in range(nAnnotation):
    annotation_path = annotation_path_list[iAnnotation]
    structure_table = structure_table_list[iAnnotation]
    structure_name = os.path.splitext(os.path.basename(structure_path_list[iAnnotation]))[0]

    print(annotation_path)

    # Compute volumes
    output_table = subjectPath2volumeTable(annotation_path)

    # Add structure names to output table by using common volume integers
    output_table = pd.merge(left=output_table,
                            right=structure_table_list[iAnnotation].loc[:, ['VolumeInteger', 'name']],
                            left_on='VolumeInteger',
                            right_on='VolumeInteger')

    # Add metadata
    output_table['atlas'] = structure_name

    # Add table to table list
    output_table_list.append(output_table)

output_table_all = pd.concat(output_table_list, ignore_index=True)
output_table_all.to_csv(os.path.join(analysis_path, 'reference_volumes.csv'))
reference_table = output_table_all.copy()
reference_cerebellum_table = reference_table[reference_table['atlas'] == 'Lobules-SUIT']
reference_cerebellum_table.to_csv(os.path.join(analysis_path, 'reference_cerebellum_volumes.csv'))
# reference_table[reference_table['']=='']



## Calculate volumes for all input files
output_table_list = list()
for iAnnotation in range(nAnnotation):
    input_path_list = input_path_list_list[iAnnotation]
    structure_table = structure_table_list[iAnnotation]
    structure_name = os.path.splitext(os.path.basename(structure_path_list[iAnnotation]))[0]

    for iInput, Input in enumerate(input_path_list):
        input_name = Input.split(os.sep)[-2]
        print('subject '+str(iInput)+': '+input_name)

        # Compute volumes
        output_table = subjectPath2volumeTable(Input)

        # Add structure names to output table by using common volume integers
        output_table = pd.merge(left=output_table,
                                right=structure_table_list[iAnnotation].loc[:, ['VolumeInteger', 'name']],
                                left_on='VolumeInteger',
                                right_on='VolumeInteger')

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

output_table_all.to_csv(os.path.join(analysis_path, 'all_volumes.csv'))
output_table_cerebellum = output_table_all[output_table_all['atlas'] == 'Lobules-SUIT']
output_table_cerebellum.to_csv(os.path.join(analysis_path, 'cerebellum_volumes.csv'))



## Go through each structure in volume table and summarize statistics per structure
output_table_pername_list = list()
for nameStruct in np.unique(np.array(output_table_all['name'].astype('category'))):

    output_table_pername = output_table_all[output_table_all['name'] == nameStruct]
    atlas_name = output_table_pername['atlas'].iloc[0]
    controlMeanVolume = np.mean(np.array([float(output_table_pername[output_table_pername['subject'] == 'control1']['VolumeNormalized']),
                                          float(output_table_pername[output_table_pername['subject'] == 'control2']['VolumeNormalized'])]))
    patientVolume = float(output_table_pername[output_table_pername['subject'] == 'patient']['VolumeNormalized'])
    fractionDifferenceVolume = (patientVolume - controlMeanVolume) / controlMeanVolume

    VolumeInteger_val = np.array(output_table_pername['VolumeInteger'])[0]

    output_table_pername_list.append(pd.DataFrame({'name': [nameStruct],
                                                  'atlas': [atlas_name],
                                                   'VolumeInteger': [VolumeInteger_val],
                                                  'fractionDifference': [fractionDifferenceVolume],
                                                  'controlMeanVolume': [controlMeanVolume],
                                                  'patientVolume': [patientVolume]}))

output_table_pername = pd.concat(output_table_pername_list, ignore_index=True)
output_table_pername = output_table_pername.sort_values(by='fractionDifference')
output_table_pername['relFracDiff'] = output_table_pername['fractionDifference']/float(output_table_pername[output_table_pername['name']=='mask']['fractionDifference'])
output_table_pername.to_csv(os.path.join(analysis_path, 'pername'+'_volumes.csv'))
output_table_pername_cerebellum = output_table_pername[output_table_pername['atlas'] == 'Lobules-SUIT']
output_table_pername_cerebellum.to_csv(os.path.join(analysis_path, 'pername_cerebellum_volumes.csv'))



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
VOIs = ['substantia nigra pars reticula',
        'substantia nigra pars compacta',
        'nucleus accumbens',
        'mask',
        'Left_CrusI',
        'Right_CrusI',
        'Vermis_CrusI',
        'Left_CrusII',
        'Right_CrusII',
        'Left_X',
        'Right_X']
for iVOI in range(len(VOIs)):
    output_table_all_plot = output_table_all[output_table_all['name'] == VOIs[iVOI]]
    ax = output_table_all_plot.plot.bar(x='subject',
                                        y='VolumeNormalized',
                                        rot=0,
                                        color=[[0.4, 0.4, 0.85],
                                               [0.4, 0.4, 0.85],
                                               [0.4, 0.4, 0.85],
                                               [0.85, 0.4, 0.4]])
    plt.ylabel('Volume Normalized')
    plt.title(VOIs[iVOI])
    ax.get_legend().remove()
    plt.show()
    plt.savefig(os.path.join(analysis_path,
                             'volume_barplot_'+VOIs[iVOI]+'.png'))
