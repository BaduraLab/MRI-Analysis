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
from functions import getChildStructures_mouse
from functions import zeroPadImage
import copy



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
lowdetail2_path_list = glob.glob(os.path.join(data_path, '*', '*_lowdetail2_manual.nii.gz')) + \
                       glob.glob(os.path.join(reference_path, '*', '*_lowdetail2.nii.gz'))
allen_image_path = os.path.join(reference_path, 'annotation_25_reoriented_flirted.nii.gz')
allen_template_image_path = os.path.join(reference_path, 'average_template_25_reoriented_flirted.nii.gz')
reference_structure_path = os.path.join(reference_path, 'structure_graph_plus.csv')
# allen_image_flirted = os.path.join(allen_path, 'annotation_25_to_AMBMC_flirted.nii.gz')
# voxel_volume = pow(0.05, 3)
# voxel_reference_volume = voxel_volume
# voxel_volume = 0.000125
# voxel_reference_volume = 1.5625e-5
include_COMcentered_images = False



# Get or define structure group id lists
allen_table = pd.read_csv(reference_structure_path)
allen_table['in_cerebellum'] = False
for iVolume in range(allen_table.shape[0]):
    if isinstance(allen_table.loc[iVolume, 'structure_id_path_custom'], str):
        allen_table.loc[iVolume, 'in_cerebellum'] = 1186 in list(map(int, allen_table.loc[iVolume, 'structure_id_path_custom'].strip('][').split(', ')))
cerebellum_ids = np.array(allen_table[allen_table['in_cerebellum']]['id_custom'])
cerebellum_ids = np.round(cerebellum_ids[~np.isnan(cerebellum_ids)]).astype(int)
sn_ids = [54, 268] # compact and reticular respectively
cerebellum_sig_ids = [821, 487] # cerebellum ids of significant structures (pValFDR_BrainNormalized)
expression_selection_ids = [502, 619, 460, 473, 54, 268, 852, 916, 75, 1241, 2001, 607] # respectively: ['Periaqueductal gray', 'Parabigeminal nucleus', 'Midbrain reticular nucleus', 'Interpeduncular nucleus', 'Substantia nigra, reticular part', 'Substantia nigra, compact part', 'Parabrachial nucleus', 'Cuneiform nucleus', 'Pedunculopontine nucleus', 'Ventral tegmental area', 'Anterior pretectal nucleus']
# group_names_list = ['ci', 'si', 'ceSig', 'expSel', 'all']

########################################################################################################################

# group_ids_list = [list(cerebellum_ids), sn_ids, cerebellum_sig_ids, expression_selection_ids, list(np.unique(allen_table['id_custom']))]




# group_names_list = ['ci', 'expSel']
# group_ids_list = [list(cerebellum_ids), expression_selection_ids]



group_names_list = ['mask']
group_ids_list = [[855]]



mb_ids = [834] # [660, 715, 504, 1317, 508, 140, 418, 476, 299, 850, 782, 371, 599]
# mb_ids = [834, 660, 715, 504, 1317, 508, 140, 418, 476, 299, 850, 782, 371, 599]

########################################################################################################################


# # Isolate structures for allen atlas [ce, sn, ce_sig, exp_sel]
# allen_image = nib.load(allen_image_path)
# allen = allen_image.get_fdata()
# allen = np.round(allen).astype(int)
# allen_template_image = nib.load(allen_template_image_path)
# allen_template = allen_template_image.get_fdata()
# for iGroup in range(len(group_names_list)):
#     group_dir = os.path.join(reference_path, group_names_list[iGroup])
#     if not os.path.exists(group_dir):
#         os.makedirs(group_dir)
#
#     child = getChildStructures_mouse(group_ids_list[iGroup], allen_table)
#     allen_in_group = np.isin(allen, child) # alter these so that they also contain all lower level structures
#
#     allen_group = allen * allen_in_group
#     allen_template_group = allen_template * allen_in_group
#
#     save_image(allen_group, allen_image, os.path.join(group_dir, 'annotation_50_reoriented_' + '_' + group_names_list[iGroup] + '.nii.gz'))
#     save_image(allen_template_group, allen_template_image, os.path.join(group_dir, 'average_template_50_reoriented_' + '_' + group_names_list[iGroup] + '.nii.gz'))
#
#     for iID in range(len(group_ids_list[iGroup])):
#         child = getChildStructures_mouse([group_ids_list[iGroup][iID]], allen_table)
#         allen_in_ID = np.isin(allen, child)
#
#         allen_ID = allen * allen_in_ID
#         allen_template_ID = allen_template * allen_in_ID
#
#         if group_names_list[iGroup] == 'all':
#             acr = allen_table.loc[allen_table['id_custom'] == group_ids_list[iGroup][iID], 'name'].iloc[0]
#             acr = acr.replace('/', '')
#         else:
#             acr = allen_table.loc[allen_table['id_custom'] == group_ids_list[iGroup][iID], 'acronym'].iloc[0]
#         save_image(allen_ID, allen_image, os.path.join(group_dir, 'annotation_50_reoriented_' + '_' + acr + '.nii.gz'))
#         save_image(allen_template_ID, allen_template_image, os.path.join(group_dir, 'average_template_50_reoriented_' + '_' + acr + '.nii.gz'))



# Isolate structures for each subject [ce, sn, ce_sig, exp_sel] for both inmasked template and annotation
# (create separate mask for each individual structure?)
input_list = glob.glob(os.path.join(data_path, '*'))
# input_list = ['Data\\Mouse\\Processed_Old\\WT_50', 'Data\\Mouse\\Processed_Old\\KO_6',
#               'Data\\Mouse\\Processed_Old\\KO_2b',
#               'Data\\Mouse\\Processed_Old\\KO_3A_2']
# input_list = input_list[0:1]
for iCase in [0,1,2,3,4]:
    for Path in input_list:
        print(Path)
        subject = Path.split(os.sep)[-1]
        print(f'subject = {subject}')

        # Change to flirted path if you want flirted annotations ################
        # annotation_path = glob.glob(os.path.join(Path, 'allen_annotation_invsynned_to_*_flirted.nii.gz'))[0]
        # template_path = glob.glob(os.path.join(Path, 'FLIRT', '*inmasked_flirted.nii.gz'))[0]
        # allen_annotation_invsynned_to_KO_3A_9_adjusted_cerebellum_lobular_flirtedRigid.nii.gz
        # KO_3A_9_inmasked_flirtedRigid.nii.gz

        ## define isolation paths for different transformation levels
        if iCase == 0:
            transform_level = 'native'
            annotation_path = glob.glob(os.path.join(Path, 'allen_annotation_invsynned_to_' + subject + '.nii.gz'))[0]
            template_path = glob.glob(os.path.join(Path, 'FLIRT', subject + '_inmasked.nii.gz'))[0]
        elif iCase == 1:
            transform_level = 'flirtedRigid'
            annotation_path = glob.glob(os.path.join(Path, 'retransform', 'allen_annotation_invsynned_to_' + subject + '_' + transform_level + '.nii.gz'))[0]
            template_path = glob.glob(os.path.join(Path, 'retransform', subject + '_' + transform_level + '.nii.gz'))[0]
        elif iCase == 2:
            transform_level = 'flirted'
            annotation_path = glob.glob(os.path.join(Path, 'retransform', 'allen_annotation_invsynned_to_' + subject + '_' + transform_level + '.nii.gz'))[0]
            template_path = glob.glob(os.path.join(Path, 'retransform', subject + '_' + transform_level + '.nii.gz'))[0]
        elif iCase == 3:
            transform_level = 'synned'
            annotation_path = glob.glob(os.path.join(Path, 'retransform', 'allen_annotation_invsynned_to_' + subject + '_' + transform_level + '.nii.gz'))[0]
            template_path = glob.glob(os.path.join(Path, 'retransform', subject + '_' + transform_level + '.nii.gz'))[0]
        else:
            transform_level = 'native-adjusted'
            annotation_path = glob.glob(os.path.join(Path, 'allen_annotation_invsynned_to_' + subject + '_adjusted_cerebellum_lobular.nii.gz'))[0]
            template_path = glob.glob(os.path.join(Path, 'FLIRT', subject + '_inmasked.nii.gz'))[0]

        Path_iso = os.path.join(Path, 'Isolations_' + transform_level)
        if not os.path.exists(Path_iso):
            os.makedirs(Path_iso)

        annotation_image = nib.load(annotation_path)
        annotation = annotation_image.get_fdata()
        annotation = np.round(annotation).astype(int)
        template_image = nib.load(template_path)
        template = template_image.get_fdata()

        for iGroup in range(len(group_names_list)):
            group_dir = os.path.join(Path_iso, group_names_list[iGroup])
            if not os.path.exists(group_dir):
                os.makedirs(group_dir)

            child = getChildStructures_mouse(group_ids_list[iGroup], allen_table)
            annotation_in_group = np.isin(annotation, child)

            # Isolate structures
            annotation_group = annotation * annotation_in_group
            template_group = template * annotation_in_group

            # Save non-COM-centered isolated structures
            save_image(annotation_group, annotation_image,
                       os.path.join(group_dir, annotation_path.split('.')[0].split(os.sep)[-1] + '_' + group_names_list[iGroup] + '.nii.gz'))
            save_image(template_group, template_image,
                       os.path.join(group_dir, template_path.split('.')[0].split(os.sep)[-1] + '_' + group_names_list[iGroup] + '.nii.gz'))

            group_xyz_maxCross_table_list = []
            for iID in range(len(group_ids_list[iGroup])):
                iID_id_custom = group_ids_list[iGroup][iID]
                iID_name = allen_table.loc[allen_table['id_custom'] == iID_id_custom, 'name'].iloc[0]
                iID_acronym = allen_table.loc[allen_table['id_custom'] == iID_id_custom, 'acronym'].iloc[0]
                child = getChildStructures_mouse(iID_id_custom, allen_table)

                annotation_in_ID = np.isin(annotation, child)
                if np.any(annotation_in_ID):
                    annotation_ID = annotation * annotation_in_ID
                    template_ID = template * annotation_in_ID

                    if group_names_list[iGroup] == 'all':
                        # acr = allen_table.loc[allen_table['id_custom'] == group_ids_list[iGroup][iID], 'name'].iloc[0]
                        # acr = acr.replace('/', '')
                        acr = allen_table.loc[allen_table['id_custom'] == group_ids_list[iGroup][iID], 'acronym'].iloc[0]
                        acr = acr.replace('/', '')
                    else:
                        acr = allen_table.loc[allen_table['id_custom'] == group_ids_list[iGroup][iID], 'acronym'].iloc[0]
                    save_image(annotation_ID, annotation_image, os.path.join(group_dir, annotation_path.split('.')[0].split(os.sep)[-1] + '_' + acr + '.nii.gz'))
                    save_image(template_ID, template_image,
                               os.path.join(group_dir, template_path.split('.')[0].split(os.sep)[-1] + '_' + acr + '.nii.gz'))

                    if include_COMcentered_images:
                        ## Only get COM and COMcenter and LR separate and save crossection table if not 'all'
                        if group_names_list[iGroup] != 'all':
                            # get COM original image
                            qform = annotation_image.get_qform()
                            voxel_volume = abs(np.linalg.det(qform))
                            COM_voxCoords = np.mean(np.array(np.where(annotation_in_ID > 0)), axis=1)
                            COM = np.matmul(qform, np.concatenate([COM_voxCoords, [1]]))

                            # get max crossection coronal, sagittal and axial coordinates (also voxCoords)
                            # Go through x-y-z, compute crossection, afterwards take maximum and remember x-y-z coordinate
                            image_shape = annotation_in_ID.shape

                            area_x_nVoxel = np.zeros([image_shape[0], 1])
                            area_x = np.zeros([image_shape[0], 1])
                            area_x_center_voxCoord_list = []
                            for x in range(image_shape[0]):
                                annotation_in_ID_slice = np.squeeze(annotation_in_ID[x, :, :])
                                area_x_nVoxel[x] = np.sum(annotation_in_ID_slice)
                                area_x[x] = area_x_nVoxel[x] * voxel_volume
                                area_x_center_voxCoord_list.append(np.mean(np.array(np.where(annotation_in_ID_slice > 0)), axis=1))
                            voxCoord_x = np.argmax(area_x)
                            area_x_center_voxCoord = np.array([voxCoord_x, area_x_center_voxCoord_list[voxCoord_x][0], area_x_center_voxCoord_list[voxCoord_x][1]])
                            area_x_center_Coord = np.matmul(qform, np.concatenate([area_x_center_voxCoord, [1]]))
                            area_x_center_Coord = area_x_center_Coord[0:3]
                            area_x_center_mmCoord = list(area_x_center_Coord / 1e3)
                            area_x_center_Coord = list(area_x_center_Coord)
                            area_x_nVoxel_max = area_x_nVoxel[voxCoord_x]
                            area_x_max = area_x[voxCoord_x]
                            template_x_minmax = [np.min(template[voxCoord_x, :, :]), np.max(template[voxCoord_x, :, :])]

                            area_y_nVoxel = np.zeros([image_shape[1], 1])
                            area_y = np.zeros([image_shape[1], 1])
                            area_y_center_voxCoord_list = []
                            for y in range(image_shape[1]):
                                annotation_in_ID_slice = np.squeeze(annotation_in_ID[:, y, :])
                                area_y_nVoxel[y] = np.sum(annotation_in_ID_slice)
                                area_y[y] = area_y_nVoxel[y] * voxel_volume
                                area_y_center_voxCoord_list.append(np.mean(np.array(np.where(annotation_in_ID_slice > 0)), axis=1))
                            voxCoord_y = np.argmax(area_y)
                            area_y_center_voxCoord = np.array([area_y_center_voxCoord_list[voxCoord_y][0], voxCoord_y, area_y_center_voxCoord_list[voxCoord_y][1]])
                            area_y_center_Coord = np.matmul(qform, np.concatenate([area_y_center_voxCoord, [1]]))
                            area_y_center_Coord = area_y_center_Coord[0:3]
                            area_y_center_mmCoord = list(area_y_center_Coord / 1e3)
                            area_y_center_Coord = list(area_y_center_Coord)
                            area_y_nVoxel_max = area_y_nVoxel[voxCoord_y][0]
                            area_y_max = area_y[voxCoord_y]
                            template_y_minmax = [np.min(template[:, voxCoord_y, :]), np.max(template[:, voxCoord_y, :])]

                            area_z_nVoxel = np.zeros([image_shape[2], 1])
                            area_z = np.zeros([image_shape[2], 1])
                            area_z_center_voxCoord_list = []
                            for z in range(image_shape[2]):
                                annotation_in_ID_slice = np.squeeze(annotation_in_ID[:, :, z])
                                area_z_nVoxel[z] = np.sum(annotation_in_ID_slice)
                                area_z[z] = area_z_nVoxel[z] * voxel_volume
                                area_z_center_voxCoord_list.append(np.mean(np.array(np.where(annotation_in_ID_slice > 0)), axis=1))
                            voxCoord_z = np.argmax(area_z)
                            area_z_center_voxCoord = np.array([area_z_center_voxCoord_list[voxCoord_z][0], area_z_center_voxCoord_list[voxCoord_z][1], voxCoord_z])
                            area_z_center_Coord = np.matmul(qform, np.concatenate([area_z_center_voxCoord, [1]]))
                            area_z_center_Coord = area_z_center_Coord[0:3]
                            area_z_center_mmCoord = list(area_z_center_Coord / 1e3)
                            area_z_center_Coord = list(area_z_center_Coord)
                            area_z_nVoxel_max = area_z_nVoxel[voxCoord_z]
                            area_z_max = area_z[voxCoord_z]
                            template_z_minmax = [np.min(template[:, :, voxCoord_z]), np.max(template[:, :, voxCoord_z])]

                            # Add table entry to list
                            group_xyz_maxCross_table_list.append(pd.DataFrame({'name': [iID_name],
                                                                              'acronym': [iID_acronym],
                                                                              'id_custom': [iID_id_custom],
                                                                              'voxCoord_x': [voxCoord_x],
                                                                              'area_x_center_voxCoord': [area_x_center_voxCoord],
                                                                              'area_x_center_Coord': [area_x_center_Coord],
                                                                              'area_x_center_mmCoord': [area_x_center_mmCoord],
                                                                              'area_x_nVoxel_max': [area_x_nVoxel_max],
                                                                              'area_x_max': [area_x_max],
                                                                              'template_x_minmax': [template_x_minmax],
                                                                              'voxCoord_y': [voxCoord_y],
                                                                              'area_y_center_voxCoord': [area_y_center_voxCoord],
                                                                              'area_y_center_Coord': [area_y_center_Coord],
                                                                              'area_y_center_mmCoord': [area_y_center_mmCoord],
                                                                              'area_y_nVoxel_max': [area_y_nVoxel_max],
                                                                              'area_y_max': [area_y_max],
                                                                              'template_y_minmax': [template_y_minmax],
                                                                              'voxCoord_z': [voxCoord_z],
                                                                              'area_z_center_voxCoord': [area_z_center_voxCoord],
                                                                              'area_z_center_Coord': [area_z_center_Coord],
                                                                              'area_z_center_mmCoord': [area_z_center_mmCoord],
                                                                              'area_z_nVoxel_max': [area_z_nVoxel_max],
                                                                              'area_z_max': [area_z_max],
                                                                              'template_z_minmax': [template_z_minmax]}))


                            # left-right and COM to [0,0,0]
                            for iLR in range(3):
                                # Get annotation_in_ID_LR
                                if iLR != 0:
                                    where_array = np.array(np.where(annotation_in_ID))
                                    annotation_in_ID_LR = annotation_in_ID.copy()
                                    where_array_coords = np.zeros(where_array.shape)
                                    for iC in range(where_array.shape[1]):
                                        where_array_coords[:, iC] = np.matmul(qform, np.concatenate([where_array[:, iC], [1]]))[0:3]
                                        if iLR == 1:
                                            if where_array_coords[0, iC] > COM[0]:
                                                annotation_in_ID_LR[tuple(where_array[:, iC])] = 0
                                        elif iLR == 2:
                                            if where_array_coords[0, iC] < COM[0]:
                                                annotation_in_ID_LR[tuple(where_array[:, iC])] = 0
                                else:
                                    annotation_in_ID_LR = 1
                                annotation_ID_zeroPadded = annotation_ID * annotation_in_ID_LR
                                template_ID_zeroPadded = template_ID * annotation_in_ID_LR

                                # crop
                                annotation_ID_zeroPadded, crop_index = zeroPadImage(annotation_ID_zeroPadded, annotation_ID, 0.2)
                                template_ID_zeroPadded, crop_index = zeroPadImage(template_ID_zeroPadded, annotation_ID, 0.2)

                                # get COM of new image and set to [0, 0, 0]
                                qform_new = annotation_image.get_qform()
                                COM_new = np.mean(np.array(np.where(annotation_ID_zeroPadded > 0)), axis=1)
                                COM_new = np.matmul(qform_new, np.concatenate([COM_new, [1]]))
                                qform_new[0, 3] = qform_new[0, 3] - COM_new[0]
                                qform_new[1, 3] = qform_new[1, 3] - COM_new[1]
                                qform_new[2, 3] = qform_new[2, 3] - COM_new[2]
                                annotation_image_COM = copy.deepcopy(annotation_image)
                                annotation_image_COM.set_qform(qform_new)
                                template_image_COM = copy.deepcopy(template_image)
                                template_image_COM.set_qform(qform_new)

                                if iLR == 0:
                                    save_image(annotation_ID_zeroPadded, annotation_image_COM, os.path.join(group_dir, annotation_path.split('.')[0].split(os.sep)[-1] + '_' + acr + '_COMcentered.nii.gz'))
                                    save_image(template_ID_zeroPadded, template_image_COM,
                                               os.path.join(group_dir, template_path.split('.')[0].split(os.sep)[-1] + '_' + acr + '_COMcentered.nii.gz'))
                                elif iLR == 1:
                                    save_image(annotation_ID_zeroPadded, annotation_image_COM, os.path.join(group_dir,
                                                                                             annotation_path.split('.')[
                                                                                                 0].split(os.sep)[
                                                                                                 -1] + '_' + acr + '_COMcentered_left.nii.gz'))
                                    save_image(template_ID_zeroPadded, template_image_COM,
                                               os.path.join(group_dir, template_path.split('.')[0].split(os.sep)[
                                                   -1] + '_' + acr + '_COMcentered_left.nii.gz'))
                                elif iLR == 2:
                                    save_image(annotation_ID_zeroPadded, annotation_image_COM, os.path.join(group_dir,
                                                                                             annotation_path.split('.')[
                                                                                                 0].split(os.sep)[
                                                                                                 -1] + '_' + acr + '_COMcentered_right.nii.gz'))
                                    save_image(template_ID_zeroPadded, template_image_COM,
                                               os.path.join(group_dir, template_path.split('.')[0].split(os.sep)[
                                                   -1] + '_' + acr + '_COMcentered_right.nii.gz'))

            if include_COMcentered_images:
                if group_names_list[iGroup] != 'all':
                    group_xyz_maxCross_table = pd.concat(group_xyz_maxCross_table_list, ignore_index=True)
                    group_xyz_maxCross_table.to_csv(os.path.join(group_dir, group_names_list[iGroup] + '.csv'))


# Manually adjusted midbrain isolation
lowdetail2_template_path_list = [os.path.join(data_path, 'WT_30_female', 'WT_30_female_reoriented.nii.gz'),
                                 os.path.join(reference_path, 'average_template_25_reoriented.nii.gz')]
lowdetail2_annotation_path_list = [glob.glob(os.path.join(data_path, 'WT_30_female', '*WT_30_female_*_lowdetail2_manual.nii.gz'))[0],
                                   os.path.join(reference_path, 'annotation_25_reoriented_lowdetail2.nii.gz')]
annotation_path_list = [glob.glob(os.path.join(data_path, 'WT_30_female', '*WT_30_female_*_lobular.nii.gz'))[0],
                        os.path.join(reference_path, 'annotation_25_reoriented.nii.gz')]
for iLowdetail2 in range(len(lowdetail2_template_path_list)):
    # Define input and output paths
    lowdetail2_template_path = lowdetail2_template_path_list[iLowdetail2]
    lowdetail2_annotation_path = lowdetail2_annotation_path_list[iLowdetail2]
    annotation_path = annotation_path_list[iLowdetail2]
    lowdetail2_template_masked_path = lowdetail2_template_path.split('.')[0] + '_mbmasked.nii.gz'
    lowdetail2_annotation_masked_path = lowdetail2_annotation_path.split('.')[0] + '_mbmasked.nii.gz'
    annotation_masked_path = annotation_path.split('.')[0] + '_mbmasked.nii.gz'

    # Load images
    print(f'Loading {lowdetail2_template_path}')
    template_image = nib.load(lowdetail2_template_path)
    template = template_image.get_fdata()

    print(f'Loading {lowdetail2_annotation_path}')
    annotation_lowdetail2_image = nib.load(lowdetail2_annotation_path)
    annotation_lowdetail2 = annotation_lowdetail2_image.get_fdata()
    annotation_lowdetail2 = np.round(annotation_lowdetail2).astype(int)

    print(f'Loading {annotation_path}')
    annotation_image = nib.load(annotation_path)
    annotation = annotation_image.get_fdata()
    annotation = np.round(annotation).astype(int)

    # Mask images with midbrain annotation
    midbrain_mask = np.isin(annotation_lowdetail2, mb_ids)
    template_masked = template * midbrain_mask
    annotation_lowdetail2_masked = annotation_lowdetail2 * midbrain_mask
    annotation_masked = annotation * midbrain_mask

    # Save midbrain masked images
    save_image(template_masked, template_image, lowdetail2_template_masked_path)
    save_image(annotation_lowdetail2_masked, annotation_lowdetail2_image, lowdetail2_annotation_masked_path)
    save_image(annotation_masked, annotation_image, annotation_masked_path)
