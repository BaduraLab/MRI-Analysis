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
                        glob.glob(os.path.join(data_path, '*', '*annotation_orsuit_thrarg_adjusted.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_orsuit_thrarg.nii.gz'))]
# glob.glob(os.path.join(data_path, '*', '*annotation_SEG.nii.gz'))
annotation_path = os.path.join(reference_path,
                                 'suit',
                                 'atlasesSUIT',
                                 'Lobules-SUIT.nii')
structure_path = os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT.csv')
reference_volumes_path = os.path.join(analysis_path, 'reference_volumes.csv')
structure_name_list = ['adjusted_lobular', 'adjusted', 'automatic']

# Follows
nAnnotation = len(input_path_list_list)
structure_table = pd.read_csv(structure_path)
reference_table = pd.read_csv(reference_volumes_path)


## For each nifti calculate volumes, then later merge each of these tables
# load unadjusted orsuit volume
# load cerebellar adjusted orsuit volume
# load final manually adjusted orsuit volume
# load SEG volume
# calculate volumes
# calculate volume differences between successive stages
# order according to highest relative volume differences from 2nd to 3rd stage

## Calculate volumes for all files
output_table_list = list()
for iAnnotation in range(nAnnotation):
    input_path_list = input_path_list_list[iAnnotation]
    structure_name = structure_name_list[iAnnotation]
    print(annotation_path)

    for iInput in range(len(input_path_list)):
        input_path = input_path_list[iInput]
        subject = input_path.split(os.sep)[-2]
        print(input_path)
        print(subject)

        # Compute volumes
        output_table = subjectPath2volumeTable(input_path)

        # Add structure names to output table by using common volume integers
        output_table = pd.merge(left=output_table,
                                right=structure_table.loc[:, ['VolumeInteger', 'name']],
                                left_on='VolumeInteger',
                                right_on='VolumeInteger')

        # Add metadata
        output_table['annotation'] = structure_name
        output_table['subject'] = subject

        # Add reference volume
        reference_table = reference_table.rename(columns={'Volume': 'rVolume'})
        output_table = pd.merge(left=output_table, right=reference_table[['name', 'VolumeInteger', 'rVolume']],
                                left_on=['name', 'VolumeInteger'], right_on=['name', 'VolumeInteger'])
        output_table['nVolume'] = output_table['Volume'] / output_table['rVolume']

        # Add table to table list
        output_table_list.append(output_table)

output_table_all = pd.concat(output_table_list, ignore_index=True)
output_table_pivoted = pd.pivot_table(output_table_all,
               index=['subject', 'VolumeInteger', 'name'],
               columns='annotation',
               values='nVolume').\
    reindex(columns=['automatic', 'adjusted', 'adjusted_lobular'])
output_table_pivoted = output_table_pivoted.sort_values(by=['VolumeInteger', 'name', 'subject'])
output_table_pivoted.to_csv(os.path.join(analysis_path, 'manual_analysis_volumes.csv'))
# reference_table = output_table_all.copy()
# reference_cerebellum_table = reference_table[reference_table['atlas'] == 'Lobules-SUIT']
# reference_cerebellum_table.to_csv(os.path.join(analysis_path, 'reference_cerebellum_volumes.csv'))
# reference_table[reference_table['']=='']



# go through each structure and count number of removed and added voxels,
# which gives an actual idea of annotation instead of