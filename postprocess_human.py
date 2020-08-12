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

# Define paths
data_path = os.path.join('Data', 'Human', 'Processed')
reference_path = os.path.join('Data', 'Human', 'Reference')
analysis_path = os.path.join('Data', 'Human', 'Analysis')
structure_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT.txt'),
                       os.path.join(reference_path, 'subcortical', 'subcortical.csv'),
                       os.path.join(reference_path, 'CerebrA', 'CerebrA.csv'),
                       os.path.join(reference_path, 'CerebrA', 'mask.csv')]



## Convert orsuit to correct orientation (SUIT method output)
orsuit_path_list = glob.glob(os.path.join(data_path, '*', '*annotation_orsuit_thrarg.nii'))
for iPath in orsuit_path_list:
    print(iPath)
    image = nib.load(iPath)
    print(nib.aff2axcodes(image.affine))
    image_reoriented = nib.as_closest_canonical(image)
    nib.save(image_reoriented, iPath+'.gz')



## Create masks
# Define paths
input_path_list_list = [glob.glob(os.path.join(data_path, '*', '*annotation_suit_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_orsuit_thrarg.nii.gz'))]

# Follows
structure_table_list = [pd.read_csv(structure_path_list[0],
                                    delim_whitespace=True,
                                    names=['VolumeInteger', 'name', 'ID']),
                        pd.read_csv(structure_path_list[1],
                                    names=['VolumeInteger', 'acronym', 'name']),
                        pd.read_csv(structure_path_list[2]),
                        pd.read_csv(structure_path_list[3])]



# Create suit mask for manual adjustment
for iPath in sum(input_path_list_list, []):
    print(iPath)

    # Get mask output path
    mask_path = iPath.split('_thrarg')[0] + 'mask_thrarg.nii.gz'
    print(mask_path)

    # Load image
    input_image = nib.load(iPath)
    input = input_image.get_fdata()

    mask = input > 0
    mask_image = nib.Nifti1Image(mask,
                                 input_image.affine,
                                 input_image.header)
    nib.save(mask_image, mask_path)
