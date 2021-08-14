import os
import glob
import nibabel as nib
import numpy as np
from functions import save_image

"""
Split KO, patient and control means as well as effect sizes into 2 files for positive and negative values.
This works better for visualization in mricrogl.
"""

# paths to files - mouse
data_mouse_path = os.path.join('Data', 'Mouse', 'Processed_Old')
reference_mouse_path = os.path.join('Data', 'Mouse', 'Reference')
analysis_mouse_path = os.path.join('Data', 'Mouse', 'Analysis')
SD_mouse_path = os.path.join(analysis_mouse_path, 'SD_JDet')
data_mouse_file_path_list = glob.glob(os.path.join(SD_mouse_path, 'SD_syn_defField_magnitude_meanWT_*.nii.gz')) + \
                            glob.glob(os.path.join(SD_mouse_path, 'SD_syn_defField_magnitude_meanKO_*.nii.gz')) + \
                            glob.glob(os.path.join(SD_mouse_path, 'SD_syn_defField_magnitude_CohenD_*.nii.gz'))
data_human_path = os.path.join('Data', 'Human', 'Processed')
reference_human_path = os.path.join('Data', 'Human', 'Reference')
analysis_human_path = os.path.join('Data', 'Human', 'Analysis')
SD_human_path = os.path.join(analysis_human_path, 'SD_JDet')
POS_mouse_path = os.path.join(SD_mouse_path, 'POS')
if not os.path.exists(POS_mouse_path):
    os.makedirs(POS_mouse_path)
NEG_mouse_path = os.path.join(SD_mouse_path, 'NEG')
if not os.path.exists(NEG_mouse_path):
    os.makedirs(NEG_mouse_path)

# paths to files - human
data_human_file_path_list = glob.glob(os.path.join(SD_human_path, 'SD_syn_defField_magnitude_meanControl_*.nii.gz')) + \
                            glob.glob(os.path.join(SD_human_path, 'SD_syn_defField_magnitude_Patient_*.nii.gz')) + \
                            glob.glob(os.path.join(SD_human_path, 'SD_syn_defField_magnitude_CohenD_*.nii.gz'))
POS_human_path = os.path.join(SD_human_path, 'POS')
if not os.path.exists(POS_human_path):
    os.makedirs(POS_human_path)
NEG_human_path = os.path.join(SD_human_path, 'NEG')
if not os.path.exists(NEG_human_path):
    os.makedirs(NEG_human_path)

# paths to files - combine human and mouse
data_path_list = data_mouse_file_path_list + data_human_file_path_list

# for each of the files, separate them into positive and negative volumes, taking the absolute of negative values
for data_path in data_path_list:
    data_image = nib.load(data_path)
    data = data_image.get_fdata()

    data_pos = data.copy()
    data_neg = data.copy()

    data_pos[data < 0] = 0
    data_neg[data > 0] = 0
    data_neg = np.abs(data_neg)

    # split between human and mouse
    if data_path.split(os.sep)[-4] == 'Human':
        data_pos_path = os.path.join(POS_human_path, data_path.split(os.sep)[-1].split('.')[0] + '_pos.nii.gz')
        save_image(data_pos, data_image, data_pos_path)
        data_neg_path = os.path.join(NEG_human_path, data_path.split(os.sep)[-1].split('.')[0] + '_neg.nii.gz')
        save_image(data_neg, data_image, data_neg_path)
    elif data_path.split(os.sep)[-4] == 'Mouse':
        data_pos_path = os.path.join(POS_mouse_path, data_path.split(os.sep)[-1].split('.')[0] + '_pos.nii.gz')
        save_image(data_pos, data_image, data_pos_path)
        data_neg_path = os.path.join(NEG_mouse_path, data_path.split(os.sep)[-1].split('.')[0] + '_neg.nii.gz')
        save_image(data_neg, data_image, data_neg_path)