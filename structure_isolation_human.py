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
from functions import zeroPadImage
import copy

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
annotation_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT.nii.gz'),
                        os.path.join(reference_path, 'subcortical', 'prob_atlas_bilateral_thrarg_0.4.nii.gz')]
template_path_list = [os.path.join(reference_path, 'suit', 'templates', 'SUIT.nii'),
                      os.path.join(reference_path, 'subcortical', 'CIT168_T1w_700um_reoriented.nii.gz')]
structure_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT_mc.csv'),
                       os.path.join(reference_path, 'subcortical', 'subcortical.csv')]
structure_list = [pd.read_csv(structure_path_list[0]),
                  pd.read_csv(structure_path_list[1])]

# Define cerebellum and substantia nigra ids
# group_ids_list = [[25, 17, 30, 4, 9, 6, 31, 8, 5, 20, 22, 23, 19, 14, 15, 16, 2, 18, 13, 10, 21, 11, 12, 24, 1, 26, 27, 28],
                  # [7, 9, 11]]
group_ids_list = [list(structure_list[0]['VolumeInteger']), [7, 9, 11, [7, 9]]]
# group_ids_list = [[25],
#                   [7, 9, 11]]
group_name_list = ['ci', 'si']

########################################################################################################################

# # Extract cerebellum and substantia nigra for references
# for i in range(len(annotation_path_list)):
#     annotation_path = annotation_path_list[i]
#     template_path = template_path_list[i]
#
#     Path_iso = os.path.join(reference_path, 'Isolations')
#     if not os.path.exists(Path_iso):
#         os.makedirs(Path_iso)
#     group_dir = os.path.join(Path_iso, structure_name_list[i])
#     if not os.path.exists(group_dir):
#         os.makedirs(group_dir)
#
#     annotation_image = nib.load(annotation_path_list[i])
#     annotation = annotation_image.get_fdata().astype(int)
#     annotation = annotation.astype(int)
#     template_image = nib.load(template_path_list[i])
#     template = template_image.get_fdata()
#
#     annotation_in_structure = np.isin(annotation, ids_list[i])
#
#     annotation_structure = annotation * annotation_in_structure
#     template_structure = template * annotation_in_structure
#
#     save_image(annotation_structure, annotation_image,
#                os.path.join(group_dir, annotation_path_list[i].split('.')[0].split(os.sep)[-1] + '_' + structure_name_list[i] + '.nii.gz'))
#     save_image(template_structure, template_image,
#                os.path.join(group_dir, template_path_list[i].split('.')[0].split(os.sep)[-1] + '_' + structure_name_list[i] + '.nii.gz'))
#
#     for iID in range(len(ids_list[i])):
#
#         annotation_in_ID = np.isin(annotation, ids_list[i][iID])
#
#         ###
#
#         if np.any(annotation_in_ID):
#
#             annotation_ID = annotation * annotation_in_ID
#             template_ID = template * annotation_in_ID
#
#             acr = structure_list[i].loc[structure_list[i]['VolumeInteger'] == ids_list[i][iID], 'name'].iloc[0]
#             save_image(annotation_ID, annotation_image, os.path.join(group_dir, annotation_path.split('.')[0].split(os.sep)[-1] + '_' + acr + '.nii.gz'))
#             save_image(template_ID, template_image, os.path.join(group_dir, template_path.split('.')[0].split(os.sep)[-1] + '_' + acr + '.nii.gz'))
#





# # For each subject create cerebellum isolated files for both inmasked template and annotation
# annotation_path_list_list = [glob.glob(os.path.join(data_path, '*', '*orsuit_thrarg*lobular_flirtRigid.nii.gz')),
#                              glob.glob(os.path.join(data_path, '*', '*subcortical_thrarg_flirtRigid.nii.gz'))]
# # template_path_list = list(set(glob.glob(os.path.join(data_path, '*', '*reoriented.nii.gz'))) - \
# #                           set(glob.glob(os.path.join(data_path, '*', '*skull_reoriented.nii.gz'))))
# for iGroup in range(len(annotation_path_list_list)):
#     annotation_path_list = annotation_path_list_list[iGroup]
#     for iSubject in range(len(annotation_path_list)):
#         annotation_path = annotation_path_list[iSubject]
#         template_path = os.path.join(os.path.split(annotation_path)[0], annotation_path.split(os.sep)[-1].split('_')[0] + '_reoriented.nii.gz')
#         # template_path = template_path_list[iSubject]
#
#
#         Path_iso = os.path.join(os.path.split(annotation_path)[0].split('.')[0], 'Isolations')
#         if not os.path.exists(Path_iso):
#             os.makedirs(Path_iso)
#         group_dir = os.path.join(Path_iso, structure_name_list[iGroup])
#         if not os.path.exists(group_dir):
#             os.makedirs(group_dir)
#
#         annotation_image = nib.load(annotation_path)
#         annotation = annotation_image.get_fdata()
#         annotation = annotation.astype(int)
#         template_image = nib.load(template_path)
#         template = template_image.get_fdata()
#
#         annotation_in_structure = np.isin(annotation, group_ids_list[iGroup])
#
#         annotation_structure = annotation * annotation_in_structure
#         template_structure = np.squeeze(template) * annotation_in_structure
#
#         save_image(annotation_structure, annotation_image,
#                    os.path.join(group_dir, annotation_path.split('.')[0].split(os.sep)[-1] + '_' + structure_name_list[iGroup] + '.nii.gz'))
#         save_image(template_structure, template_image,
#                    os.path.join(group_dir, template_path.split('.')[0].split(os.sep)[-1] + '_' + structure_name_list[iGroup] + '.nii.gz'))
#
#         group_xyz_maxCross_table_list = []
#         for iID in range(len(group_ids_list[iGroup])):
#             iID_id_custom = group_ids_list[iGroup][iID]
#             iID_name = structure_list[iGroup].loc[structure_list[iGroup]['VolumeInteger'] == iID_id_custom, 'name'].iloc[0]
#             iID_acronym = structure_list[iGroup].loc[structure_list[iGroup]['VolumeInteger'] == iID_id_custom, 'acronym'].iloc[0]
#
#             annotation_in_ID = np.isin(annotation, group_ids_list[iGroup][iID])
#
#             annotation_ID = annotation * annotation_in_ID
#             template_ID = np.squeeze(template) * annotation_in_ID
#
#             acr = structure_list[iGroup].loc[structure_list[iGroup]['VolumeInteger'] == group_ids_list[iGroup][iID], 'name'].iloc[0]
#             save_image(annotation_ID, annotation_image, os.path.join(group_dir, annotation_path.split('.')[0].split(os.sep)[-1] + '_' + acr + '.nii.gz'))
#             save_image(template_ID, template_image, os.path.join(group_dir, template_path.split('.')[0].split(os.sep)[-1] + '_' + acr + '.nii.gz'))
#
#             # get COM original image
#             qform = annotation_image.get_qform()
#             voxel_volume = abs(np.linalg.det(qform))
#             COM = np.mean(np.array(np.where(annotation_in_ID > 0)), axis=1)
#             COM = np.matmul(qform, np.concatenate([COM, [1]]))
#
#
#
#             # get max crossection coronal, sagittal and axial coordinates (also voxCoords)
#             # Go through x-y-z, compute crossection, afterwards take maximum and remember x-y-z coordinate
#             image_shape = annotation_in_ID.shape
#
#             area_x_nVoxel = np.zeros([image_shape[0], 1])
#             area_x = np.zeros([image_shape[0], 1])
#             area_x_center_voxCoord_list = []
#             for x in range(image_shape[0]):
#                 annotation_in_ID_slice = np.squeeze(annotation_in_ID[x, :, :])
#                 area_x_nVoxel[x] = np.sum(annotation_in_ID_slice)
#                 area_x[x] = area_x_nVoxel[x] * voxel_volume
#                 area_x_center_voxCoord_list.append(np.mean(np.array(np.where(annotation_in_ID_slice > 0)), axis=1))
#             voxCoord_x = np.argmax(area_x)
#             area_x_center_voxCoord = np.array(
#                 [voxCoord_x, area_x_center_voxCoord_list[voxCoord_x][0], area_x_center_voxCoord_list[voxCoord_x][1]])
#             area_x_center_Coord = np.matmul(qform, np.concatenate([area_x_center_voxCoord, [1]]))
#             area_x_center_Coord = area_x_center_Coord[0:3]
#             area_x_center_mmCoord = list(area_x_center_Coord / 1e3)
#             area_x_center_Coord = list(area_x_center_Coord)
#             area_x_nVoxel_max = area_x_nVoxel[voxCoord_x]
#             area_x_max = area_x[voxCoord_x]
#             template_x_minmax = [np.min(template[voxCoord_x, :, :]), np.max(template[voxCoord_x, :, :])]
#
#             area_y_nVoxel = np.zeros([image_shape[1], 1])
#             area_y = np.zeros([image_shape[1], 1])
#             area_y_center_voxCoord_list = []
#             for y in range(image_shape[1]):
#                 annotation_in_ID_slice = np.squeeze(annotation_in_ID[:, y, :])
#                 area_y_nVoxel[y] = np.sum(annotation_in_ID_slice)
#                 area_y[y] = area_y_nVoxel[y] * voxel_volume
#                 area_y_center_voxCoord_list.append(np.mean(np.array(np.where(annotation_in_ID_slice > 0)), axis=1))
#             voxCoord_y = np.argmax(area_y)
#             area_y_center_voxCoord = np.array(
#                 [area_y_center_voxCoord_list[voxCoord_y][0], voxCoord_y, area_y_center_voxCoord_list[voxCoord_y][1]])
#             area_y_center_Coord = np.matmul(qform, np.concatenate([area_y_center_voxCoord, [1]]))
#             area_y_center_Coord = area_y_center_Coord[0:3]



# For each subject create cerebellum isolated files for both inmasked template and annotation
annotation_path_list_list = [glob.glob(os.path.join(data_path, '*', '*orsuit_thrarg*lobular_flirtedRigid.nii.gz')),###
                             glob.glob(os.path.join(data_path, '*', '*subcortical_thrarg_flirtedRigid.nii.gz'))]###
# template_path_list = list(set(glob.glob(os.path.join(data_path, '*', '*reoriented.nii.gz'))) - \
#                           set(glob.glob(os.path.join(data_path, '*', '*skull_reoriented.nii.gz'))))
for iGroup in range(len(annotation_path_list_list)):
    annotation_path_list = annotation_path_list_list[iGroup]
    for iSubject in range(len(annotation_path_list)):
        annotation_path = annotation_path_list[iSubject]
        template_path = os.path.join(os.path.split(annotation_path)[0], annotation_path.split(os.sep)[-1].split('_')[0] + '_flirtedRigid.nii.gz')###
        # template_path = template_path_list[iSubject]


        Path_iso = os.path.join(os.path.split(annotation_path)[0].split('.')[0], 'Isolations')
        if not os.path.exists(Path_iso):
            os.makedirs(Path_iso)
        group_dir = os.path.join(Path_iso, group_name_list[iGroup])
        if not os.path.exists(group_dir):
            os.makedirs(group_dir)

        annotation_image = nib.load(annotation_path)
        annotation = annotation_image.get_fdata()
        annotation = annotation.astype(int)
        template_image = nib.load(template_path)
        template = template_image.get_fdata()

        annotation_in_structure = np.isin(annotation, group_ids_list[iGroup])

        annotation_structure = annotation * annotation_in_structure
        template_structure = np.squeeze(template) * annotation_in_structure

        save_image(annotation_structure, annotation_image,
                   os.path.join(group_dir, annotation_path.split('.')[0].split(os.sep)[-1] + '_' + group_name_list[iGroup] + '.nii.gz'))
        save_image(template_structure, template_image,
                   os.path.join(group_dir, template_path.split('.')[0].split(os.sep)[-1] + '_' + group_name_list[iGroup] + '.nii.gz'))

        group_xyz_maxCross_table_list = []
        for iID in range(len(group_ids_list[iGroup])):
            iID_id_custom = group_ids_list[iGroup][iID]

            if isinstance(iID_id_custom, list):
                iID_name = 'SN'
            else:
                iID_name = structure_list[iGroup].loc[structure_list[iGroup]['VolumeInteger'] == iID_id_custom, 'name'].iloc[0]
            # iID_acronym = structure_list[iGroup].loc[structure_list[iGroup]['VolumeInteger'] == iID_id_custom, 'acronym'].iloc[0]

            annotation_in_ID = np.isin(annotation, iID_id_custom)

            annotation_ID = annotation * annotation_in_ID
            template_ID = np.squeeze(template) * annotation_in_ID

            save_image(annotation_ID, annotation_image, os.path.join(group_dir, annotation_path.split('.')[0].split(os.sep)[-1] + '_' + iID_name + '.nii.gz'))
            save_image(template_ID, template_image, os.path.join(group_dir, template_path.split('.')[0].split(os.sep)[-1] + '_' + iID_name + '.nii.gz'))

            # get COM original image
            qform = annotation_image.get_qform()
            voxel_volume = abs(np.linalg.det(qform))
            COM = np.mean(np.array(np.where(annotation_in_ID > 0)), axis=1)
            COM = np.matmul(qform, np.concatenate([COM, [1]]))



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
            area_x_center_voxCoord = np.array(
                [voxCoord_x, area_x_center_voxCoord_list[voxCoord_x][0], area_x_center_voxCoord_list[voxCoord_x][1]])
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
            area_y_center_voxCoord = np.array(
                [area_y_center_voxCoord_list[voxCoord_y][0], voxCoord_y, area_y_center_voxCoord_list[voxCoord_y][1]])
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
            area_z_center_voxCoord = np.array(
                [area_z_center_voxCoord_list[voxCoord_z][0], area_z_center_voxCoord_list[voxCoord_z][1], voxCoord_z])
            area_z_center_Coord = np.matmul(qform, np.concatenate([area_z_center_voxCoord, [1]]))
            area_z_center_Coord = area_z_center_Coord[0:3]
            area_z_center_mmCoord = list(area_z_center_Coord / 1e3)
            area_z_center_Coord = list(area_z_center_Coord)
            area_z_nVoxel_max = area_z_nVoxel[voxCoord_z]
            area_z_max = area_z[voxCoord_z]
            template_z_minmax = [np.min(template[:, :, voxCoord_z]), np.max(template[:, :, voxCoord_z])]

            # Add table entry to list
            group_xyz_maxCross_table_list.append(pd.DataFrame({'name': [iID_name],
                                                               'VolumeInteger': [iID_id_custom],
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
                try:
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

                    print(f'iLR={iLR}')
                    if iLR == 0:
                        save_image(annotation_ID_zeroPadded, annotation_image_COM, os.path.join(group_dir,
                                                                                                annotation_path.split('.')[
                                                                                                    0].split(os.sep)[
                                                                                                    -1] + '_' + iID_name + '_COMcentered.nii.gz'))
                        save_image(template_ID_zeroPadded, template_image_COM,
                                   os.path.join(group_dir, template_path.split('.')[0].split(os.sep)[
                                       -1] + '_' + iID_name + '_COMcentered.nii.gz'))
                    elif iLR == 1:
                        save_image(annotation_ID_zeroPadded, annotation_image_COM, os.path.join(group_dir,
                                                                                                annotation_path.split('.')[
                                                                                                    0].split(os.sep)[
                                                                                                    -1] + '_' + iID_name + '_COMcentered_left.nii.gz'))
                        save_image(template_ID_zeroPadded, template_image_COM,
                                   os.path.join(group_dir, template_path.split('.')[0].split(os.sep)[
                                       -1] + '_' + iID_name + '_COMcentered_left.nii.gz'))
                    elif iLR == 2:
                        save_image(annotation_ID_zeroPadded, annotation_image_COM, os.path.join(group_dir,
                                                                                                annotation_path.split('.')[
                                                                                                    0].split(os.sep)[
                                                                                                    -1] + '_' + iID_name + '_COMcentered_right.nii.gz'))
                        save_image(template_ID_zeroPadded, template_image_COM,
                                   os.path.join(group_dir, template_path.split('.')[0].split(os.sep)[
                                       -1] + '_' + iID_name + '_COMcentered_right.nii.gz'))
                except:
                    print('exception has occurred')

            group_xyz_maxCross_table = pd.concat(group_xyz_maxCross_table_list, ignore_index=True)
            group_xyz_maxCross_table.to_csv(os.path.join(group_dir, group_name_list[iGroup] + '.csv'))


