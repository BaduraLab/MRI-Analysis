import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import csv
from scipy.stats import ttest_ind
from pathlib import Path
import numpy_indexed as npi


## Function to compute volumes for image
def subjectPath2volumeTable(subject_path):
    # Compute voxel numbers and volumes and output to table

    # Load image
    subject_image = nib.load(subject_path)
    subject = subject_image.get_fdata()

    # Get voxel volume
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


# Define
data_path = os.path.join('Data', 'Human', 'Processed')
reference_path = os.path.join('Data', 'Human', 'Reference')
analysis_path = os.path.join('Data', 'Human', 'Analysis')
input_path_list_list = [glob.glob(os.path.join(data_path, '*', '*Cerebellum*annotation*_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*prob*annotation*_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*mni*annotation*_thrarg.nii.gz'))]
structure_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT.txt'),
                       os.path.join(reference_path, 'CIT168_Reinf_Learn_v1.1.0', 'subcortical.csv'),
                       os.path.join(reference_path, 'CerebrA', 'CerebrA.csv')]
# voxel_volume = 0.000125
# annotation_path = os.path.join(reference_path, 'annotation_50_reoriented.nii.gz')

# Follows
nAnnotation = len(input_path_list_list)
structure_table_list = [pd.read_csv(structure_path_list[0],
                                    delim_whitespace=True,
                                    names=['VolumeInteger', 'acronym', 'name']),
                        pd.read_csv(structure_path_list[1]),
                        pd.read_csv(structure_path_list[2])]




output_table_list = list()
for iAnnotation in range(nAnnotation):
    input_path_list = input_path_list_list[iAnnotation]
    structure_table = structure_table_list[iAnnotation]
    structure_name = os.path.splitext(os.path.basename(structure_path_list[iAnnotation]))

    for iInput, Input in enumerate(input_path_list):
        input_name = Input.split(os.sep)[-2]
        print('subject '+str(iInput)+': '+input_name)

        # Compute volumes
        output_table = subjectPath2volumeTable(Input)

        # Add structure names to output table by using common volume integers
        output_table = pd.merge(left=output_table,
                                right=structure_table_list[iInput].loc[:, ['VolumeInteger', 'name']],
                                left_on='VolumeInteger',
                                right_on='VolumeInteger')

        # Add metadata
        output_table['subject'] = input_name
        output_table['atlas'] = structure_name

        # Add table to table list
        output_table_list.append(output_table)

output_table_all = pd.concat(output_table_list, ignore_index=True)
output_table_all.to_csv(os.path.join(analysis_path, 'all_volumes.csv'))
