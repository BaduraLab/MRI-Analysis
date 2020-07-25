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



# Define
# mouse_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse'
# allen_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse/allen'
# data_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse/Processed_New'
# analysis_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse/Analysis'
# allen_image = os.path.join(allen_path, 'annotation_25_reoriented.nii.gz')
# allen_image_flirted = os.path.join(allen_path, 'annotation_25_to_AMBMC_flirted.nii.gz')
data_path = os.path.join('Data', 'Mouse', 'Processed')
glob.glob(os.path.join(data_path, '*'))
mouse_path = '/home/enzo/Desktop/Data/Mouse'
allen_path = '/usr/local/fsl/data/standard/allen'
data_path = '/home/enzo/Desktop/Data/Mouse/Processed_New'
analysis_path = '/home/enzo/Desktop/Data/Mouse/Analysis'
allen_image = os.path.join(allen_path, 'annotation_25_reoriented.nii.gz')
allen_image_flirted = os.path.join(allen_path, 'annotation_25_to_AMBMC_flirted.nii.gz')
picture_path = '/home/enzo/Desktop'
# voxel_volume = pow(0.05, 3)
voxel_volume = 0.000125
voxel_reference_volume = 1.5625e-5

# mouse_list = os.listdir(data_path)
# nMouse = len(mouse_list)



allen_table = pd.read_csv(os.path.join(analysis_path, 'allen_table.csv'))
allen_table['in_cerebellum'] = False
for iVolume in range(allen_table.shape[0]):
    if isinstance(allen_table.loc[iVolume, 'structure_id_path'], str):
        allen_table.loc[iVolume, 'in_cerebellum'] = 512 in list(map(int, allen_table.loc[iVolume, 'structure_id_path'].strip('][').split(', ')))

cerebellum_ids = allen_table[allen_table['in_cerebellum']]['VolumeInteger']
cerebellum_ids = cerebellum_ids[~np.isnan(cerebellum_ids)].astype(int)



# Extract cerebellum for allen atlas
allen_image_nib = nib.load(allen_image_flirted)
allen_image_array = allen_image_nib.get_fdata().astype(int)

allen_image_array_in_cerebellum = np.isin(allen_image_array, cerebellum_ids)

cerebellum_voxel_number = np.sum(allen_image_array_in_cerebellum)
cerebellum_volume = cerebellum_voxel_number * voxel_reference_volume

allen_image_array_new = allen_image_array * allen_image_array_in_cerebellum

allen_image_nib_new = nib.Nifti1Image(allen_image_array_new, allen_image_nib.affine, allen_image_nib.header)
nib.save(allen_image_nib_new, os.path.join(analysis_path, 'annotation_25_reoriented_in_cerebellum.nii.gz'))



## Get list of invwarped annotation files and compute volumes for them, compiling eventually into single table
mouse_invwarped_list = glob.glob(data_path+'/*/FLIRT/WT_50_flirted.nii.gz')
mouse_table_invwarped_list = mouse_invwarped_list.copy()


warped_image_nib = nib.load(mouse_table_invwarped_list[0])
warped_image_array = warped_image_nib.get_data()

warped_image_array_new = warped_image_array * allen_image_array_in_cerebellum

warped_image_nib_new = nib.Nifti1Image(warped_image_array_new, warped_image_nib.affine, warped_image_nib.header)
nib.save(warped_image_nib_new, os.path.join(analysis_path, 'WT_50_warped_to_allen_in_cerebellum.nii.gz'))



