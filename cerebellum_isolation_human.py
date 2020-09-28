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
data_path = os.path.join('Data', 'Human', 'Processed')
glob.glob(os.path.join(data_path, '*'))
reference_path = os.path.join('Data', 'Human', 'Reference')
analysis_path = os.path.join('Data', 'Human', 'Analysis')
annotation_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT_mc.nii.gz'),
                        os.path.join(reference_path, 'subcortical', 'prob_atlas_bilateral_thrarg_0.4.nii.gz')]
template_path_list = [os.path.join(reference_path, 'suit', 'templates', 'SUIT.nii'),
                      os.path.join(reference_path, 'subcortical', 'CIT168_T1w_700um_reoriented.nii.gz')]
structure_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT_mc.csv'),
                       os.path.join(reference_path, 'subcortical', 'subcortical.csv')]
structure_list = [pd.read_csv(structure_path_list[0]),
                  pd.read_csv(structure_path_list[1])]

# Define cerebellum and substantia nigra ids
ids_list = [[25,17,30,4,9,6,31,8,5,20,22,23,19,14,15,16,2,18,13,10,21,11,12,24,1,26,27,28],
            [7,9]]
structure_name_list = ['ci','si']



# Extract cerebellum and substantia nigra for references
for i in range(len(annotation_path_list)):
    annotation_image = nib.load(annotation_path_list[i])
    annotation = annotation_image.get_fdata().astype(int)
    annotation = annotation.astype(int)
    template_image = nib.load(template_path_list[i])
    template = template_image.get_fdata()

    annotation_in_structure = np.isin(annotation, ids_list[i])

    annotation_structure = annotation * annotation_in_structure
    template_structure = template * annotation_in_structure

    save_image(annotation_structure, annotation_image, annotation_path_list[i].split('.')[0] + '_' + structure_name_list[i] + '.nii.gz')
    save_image(template_structure, template_image, template_path_list[i].split('.')[0] + '_' + structure_name_list[i] + '.nii.gz')



# For each subject create cerebellum isolated files for both inmasked template and annotation
annotation_path_list_list = [glob.glob(os.path.join(data_path, '*', '*orsuit_thrarg*lobular_mc.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*subcortical_thrarg.nii.gz'))]
template_path_list = glob.glob(os.path.join(data_path, '*', '*reoriented.nii.gz'))
for i in range(len(annotation_path_list_list)):
    annotation_path_list = annotation_path_list_list[i]
    for iSubject in range(len(annotation_path_list)):
        annotation_path = annotation_path_list[iSubject]
        template_path = template_path_list[iSubject]

        annotation_image = nib.load(annotation_path)
        annotation = annotation_image.get_fdata()
        annotation = annotation.astype(int)
        template_image = nib.load(template_path)
        template = template_image.get_fdata()

        annotation_in_structure = np.isin(annotation, ids_list[i])

        annotation_structure = annotation * annotation_in_structure
        template_structure = np.squeeze(template) * annotation_in_structure

        save_image(annotation_structure, annotation_image, annotation_path.split('.')[0]+'_'+structure_name_list[i]+'.nii.gz')
        save_image(template_structure, template_image, template_path.split('.')[0]+'_'+structure_name_list[i]+'.nii.gz')
