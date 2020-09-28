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
from functions import save_image



# Define
# mouse_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse'
# allen_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse/allen'
# data_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse/Processed_New'
# analysis_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse/Analysis'
# allen_image = os.path.join(allen_path, 'annotation_25_reoriented.nii.gz')
# allen_image_flirted = os.path.join(allen_path, 'annotation_25_to_AMBMC_flirted.nii.gz')
data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
glob.glob(os.path.join(data_path, '*'))
reference_path = os.path.join('Data', 'Mouse', 'Reference')
analysis_path = os.path.join('Data', 'Mouse', 'Analysis')
allen_image_path = os.path.join(reference_path, 'annotation_50_reoriented_mc.nii.gz')
allen_template_image_path = os.path.join(reference_path, 'average_template_50_reoriented.nii.gz')
reference_structure_path = os.path.join(reference_path, 'structure_graph_mc.csv')
# allen_image_flirted = os.path.join(allen_path, 'annotation_25_to_AMBMC_flirted.nii.gz')
voxel_volume = pow(0.05, 3)
voxel_reference_volume = voxel_volume
# voxel_volume = 0.000125
# voxel_reference_volume = 1.5625e-5



# Get cerebellum ids
allen_table = pd.read_csv(reference_structure_path)
allen_table['in_cerebellum'] = False
for iVolume in range(allen_table.shape[0]):
    if isinstance(allen_table.loc[iVolume, 'structure_id_path'], str):
        allen_table.loc[iVolume, 'in_cerebellum'] = 512 in list(map(int, allen_table.loc[iVolume, 'structure_id_path'].strip('][').split(', ')))
cerebellum_ids = allen_table[allen_table['in_cerebellum']]['id_mc']
cerebellum_ids = np.round(cerebellum_ids[~np.isnan(cerebellum_ids)]).astype(int)

# allen_table.to_csv('allen_in_cerebellum.csv')

sn_ids = [54, 268] # compact and reticular respectively



# Extract cerebellum for allen atlas
allen_image = nib.load(allen_image_path)
allen = allen_image.get_fdata()
allen = np.round(allen).astype(int)
allen_template_image = nib.load(allen_template_image_path)
allen_template = allen_template_image.get_fdata()

allen_in_cerebellum = np.isin(allen, cerebellum_ids)

# cerebellum_voxel_number = np.sum(allen_in_cerebellum)
# cerebellum_volume = cerebellum_voxel_number * voxel_reference_volume

allen_cerebellum = allen * allen_in_cerebellum
allen_template_cerebellum = allen_template * allen_in_cerebellum

save_image(allen_cerebellum, allen_image, os.path.join(reference_path, 'annotation_50_reoriented_mc_ci.nii.gz'))
save_image(allen_template_cerebellum, allen_template_image, os.path.join(reference_path, 'average_template_50_reoriented_ci.nii.gz'))

allen_in_sn = np.isin(allen, sn_ids)

allen_sn = allen * allen_in_sn
allen_template_sn = allen_template * allen_in_sn

save_image(allen_sn, allen_image, os.path.join(reference_path, 'annotation_50_reoriented_mc_si.nii.gz'))
save_image(allen_template_sn, allen_template_image, os.path.join(reference_path, 'average_template_50_reoriented_si.nii.gz'))



# For each subject create cerebellum isolated files for both inmasked template and annotation
input_list = glob.glob(os.path.join(data_path, '*'))
# input_list = ['Data\\Mouse\\Processed_Old\\WT_50', 'Data\\Mouse\\Processed_Old\\KO_6',
#               'Data\\Mouse\\Processed_Old\\KO_2b',
#               'Data\\Mouse\\Processed_Old\\KO_3A_2']
for Path in input_list:
    print(Path)

    annotation_path = glob.glob(os.path.join(Path, '*invsynned*cerebellum_lobular_mc.nii.gz'))[0]
    template_path = glob.glob(os.path.join(Path, 'FLIRT', '*inmasked.nii.gz'))[0]

    annotation_image = nib.load(annotation_path)
    annotation = annotation_image.get_fdata()
    annotation = np.round(annotation).astype(int)
    template_image = nib.load(template_path)
    template = template_image.get_fdata()

    annotation_in_cerebellum = np.isin(annotation, cerebellum_ids)

    annotation_ci = annotation * annotation_in_cerebellum
    template_ci = template * annotation_in_cerebellum

    save_image(annotation_ci, annotation_image, annotation_path.split('.')[0]+'_ci.nii.gz')
    save_image(template_ci, template_image, template_path.split('.')[0]+'_ci.nii.gz')



    annotation_in_sn = np.isin(annotation, sn_ids)

    annotation_si = annotation * annotation_in_sn
    template_si = template * annotation_in_sn

    save_image(annotation_si, annotation_image, annotation_path.split('.')[0]+'_si.nii.gz')
    save_image(template_si, template_image, template_path.split('.')[0]+'_si.nii.gz')
