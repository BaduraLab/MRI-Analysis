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



## Define
# Function to compute volumes for image
def subjectPath2volumeTable(subject_path):
    # Compute voxel numbers and volumes and output to table

    # Load image
    subject_image = nib.load(subject_path)
    subject = subject_image.get_fdata()

    # Get voxel volume
    print(subject_image.header['pixdim'][1:4])
    voxel_volume = np.prod(subject_image.header['pixdim'][1:4]) # should be in mm^3
    print('voxel volume = '+str(voxel_volume)+'mm^3')

    # Calculate volumes
    [iVoxel, nVoxel] = np.unique(np.int64(np.round(subject)),
                                 return_counts=True)
    vVoxel = nVoxel * voxel_volume
    print('total volume = '+str(np.sum(vVoxel))+'mm^3')

    # Output to DataFrame
    volume_table = pd.DataFrame(
        {'VolumeInteger': iVoxel,
         'VoxelNumber': nVoxel,
         'Volume': vVoxel})

    return volume_table

# Define paths
data_path = os.path.join('Data', 'Human', 'Processed')
reference_path = os.path.join('Data', 'Human', 'Reference')
analysis_path = os.path.join('Data', 'Human', 'Analysis')
input_path_list_list = [glob.glob(os.path.join(data_path, '*', '*annotation_orsuit_thrarg_adjusted.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_subcortical_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_CerebrA_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_mask_thrarg.nii.gz'))]
structure_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT.txt'),
                       os.path.join(reference_path, 'subcortical', 'subcortical.csv'),
                       os.path.join(reference_path, 'CerebrA', 'CerebrA.csv'),
                       os.path.join(reference_path, 'CerebrA', 'mask.csv')]

# Follows
nAnnotation = len(input_path_list_list)
structure_table_list = [pd.read_csv(structure_path_list[0],
                                    delim_whitespace=True,
                                    names=['VolumeInteger', 'name', 'ID']),
                        pd.read_csv(structure_path_list[1],
                                    names=['VolumeInteger', 'acronym', 'name']),
                        pd.read_csv(structure_path_list[2]),
                        pd.read_csv(structure_path_list[3])]



## Calculate volumes for all atlases and input files
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
output_table_all.to_csv(os.path.join(analysis_path, 'all_volumes.csv'))



## Go through each structure in volume table and summarize statistics per structure
output_table_pername_list = list()
for nameStruct in np.unique(np.array(output_table_all['name'].astype('category'))):

    output_table_pername = output_table_all[output_table_all['name'] == nameStruct]
    atlas_name = output_table_pername['atlas'].iloc[0]
    controlMeanVolume = np.mean(np.array([float(output_table_pername[output_table_pername['subject'] == 'control1']['Volume']),
                                          float(output_table_pername[output_table_pername['subject'] == 'control2']['Volume'])]))
    patientVolume = float(output_table_pername[output_table_pername['subject'] == 'patient']['Volume'])
    fractionDifferenceVolume = (patientVolume - controlMeanVolume) / controlMeanVolume

    output_table_pername_list.append(pd.DataFrame({'name': [nameStruct],
                                                  'atlas': [atlas_name],
                                                  'fractionDifference': [fractionDifferenceVolume],
                                                  'controlMeanVolume': [controlMeanVolume],
                                                  'patientVolume': [patientVolume]}))

output_table_pername = pd.concat(output_table_pername_list, ignore_index=True)
output_table_pername = output_table_pername.sort_values(by='fractionDifference')
output_table_pername['relFracDiff'] = output_table_pername['fractionDifference']/float(output_table_pername[output_table_pername['name']=='mask']['fractionDifference'])
output_table_pername.to_csv(os.path.join(analysis_path, 'pername'+'_volumes.csv'))




## Plotting
# cmap = matplotlib.cm.get_cmap('lines')
VOIs = ['substantia nigra pars reticula',
        'substantia nigra pars compacta',
        'nucleus accumbens',
        'mask']
for iVOI in range(len(VOIs)):
    output_table_all_plot = output_table_all[output_table_all['name'] == VOIs[iVOI]]
    ax = output_table_all_plot.plot.bar(x='subject',
                                        y='Volume',
                                        rot=0,
                                        color=[[0.4, 0.4, 0.85],
                                               [0.4, 0.4, 0.85],
                                               [0.85, 0.4, 0.4]])
    plt.ylabel('Volume $mm^3$')
    plt.title(VOIs[iVOI])
    ax.get_legend().remove()
    plt.show()
    plt.savefig(os.path.join(analysis_path,
                             'volume_barplot_'+VOIs[iVOI]+'.png'))
