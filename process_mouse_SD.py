import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
# import fsl.wrappers
import csv
from scipy.stats import ttest_ind
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.imwarp import DiffeomorphicMap
from dipy.align.metrics import CCMetric
import datetime
# import SimpleITK as sitk
from compress_pickle import dump, load
from functions import zeroPadImage
from functions import save_image
from functions import imageFLIRT2defField
from functions import round_half_up

# Define functions
def isolateCropStructure(annotation, qform, children, structure_input, output_path):
    # Isolate structure in reference space
    structure_template_mask = np.isin(annotation, children)
    structure_template = structure_input * structure_template_mask  # check with fsleyes

    # Crop isolated structure
    structure_template_cropped, crop_index = zeroPadImage(structure_template, structure_template, 0.1)

    # Translate cropped window for correct affine
    qform_translated = qform.copy()
    M = qform[:3, :3]
    abc = qform[:3, 3]
    crop_index_vec = M.dot(crop_index) + abc
    qform_translated[:3, 3] = crop_index_vec

    # Save cropped and translated image
    structure_template_cropped_image = nib.Nifti1Image(structure_template_cropped, qform_translated)
    nib.save(structure_template_cropped_image, output_path)
    print(output_path)

    return structure_template_cropped_image, structure_template_mask, crop_index

# Define
data_path = os.path.join('Data', 'Mouse', 'Processed')
mouse_path_list = glob.glob(os.path.join(data_path, '*'))
reference_path = os.path.join('Data', 'Mouse', 'Reference')
# average_template_50_to_AMBMC_flirted.nii.gz
# reference_template_path = os.path.join(reference_path, 'average_template_50_reoriented.nii.gz')
# reference_annotation_path = os.path.join(reference_path, 'annotation_50_reoriented.nii.gz')
reference_template_path = os.path.join(reference_path, 'average_template_50_reoriented_flirted_cropped.nii.gz')
reference_annotation_path = os.path.join(reference_path, 'annotation_50_reoriented_flirted_cropped.nii.gz')
structure_path = os.path.join(reference_path, 'structure_graph_mc.csv')
structure = pd.read_csv(structure_path)


# reference_table['in_cerebellum'] = False
# for iVolume in range(reference_table.shape[0]):
#     if isinstance(reference_table.loc[iVolume, 'structure_id_path'], str):
#         reference_table.loc[iVolume, 'in_cerebellum'] = 512 in list(map(int, reference_table.loc[iVolume, 'structure_id_path'].strip('][').split(', ')))
# reference_cerebellum_table = reference_table[reference_table['in_cerebellum']][['name', 'acronym', 'id', 'structure_id_path', 'id_custom', 'structure_id_path_custom', 'VoxelNumber', 'Volume']]
# reference_cerebellum_table.to_csv(os.path.join(analysis_path, 'reference_volumes_cerebellum.csv'))

# Get list of structure present in annotation, for each structure that is present isolate it and whatever path is below it
reference_annotation_image = nib.load(reference_annotation_path)
reference_annotation = reference_annotation_image.get_fdata()
reference_annotation = np.int64(np.round(reference_annotation))
reference_template_image = nib.load(reference_template_path)
reference_template = reference_template_image.get_fdata()

# Calculate volumes and use these to get non-zero volume children of non-zero volume structures
[reference_VolumeInteger, reference_VoxelNumber] = np.unique(reference_annotation,
                                                             return_counts=True)
reference_table = pd.DataFrame({'id_custom': reference_VolumeInteger, 'VoxelNumber': reference_VoxelNumber})
reference_table = pd.merge(left = reference_table, right=structure[['id_custom', 'name', 'structure_id_path_custom']],
                           left_on='id_custom', right_on='id_custom')
structure_id_custom_children_list = list()
for iVolume in range(reference_table.shape[0]):
    id_custom_iVolume = reference_table['id_custom'][iVolume]

    # Get all children
    structure_id_custom_children = list()
    for iVolume2 in range(reference_table.shape[0]):
        # check if current volume integer is contained as parent to find children
        structure_id_custom_path = np.array(
            list(map(int, reference_table.loc[iVolume2, 'structure_id_path_custom'].
                     strip('][').split(', '))))
        if np.any(np.isin(structure_id_custom_path, id_custom_iVolume)):
            id_custom_iVolume2 = reference_table['id_custom'][iVolume2]
            structure_id_custom_children.append(id_custom_iVolume2)

    # start_path_index = np.where(structure_id_custom_path == reference_table['id_custom'][iVolume])[0][0]
    # structure_id_custom_lower_path = structure_id_custom_path[start_path_index:]
    structure_id_custom_children_list.append(structure_id_custom_children)
reference_table['id_custom_children'] = structure_id_custom_children_list




# Take subset of reference table for testing
reference_table = reference_table[np.isin(reference_table['name'], ['Cerebellum',
                                                                    'Substantia nigra, reticular part',
                                                                    'Substantia nigra, compact part'])]
nStructure = reference_table.shape[0]



# Load mouse volumes beforehand to save reading time
mouse_image_list = list()
mouse_list = list()
mouse_annotation_image_list = list()
mouse_annotation_list = list()
mouse_defField_magnitude_5D_list = list()
for iMousePath, MousePath in enumerate(mouse_path_list):
    # Define paths
    mouse_string = MousePath.split(os.sep)[-1]
    mouse_path = os.path.join(MousePath, mouse_string + '_reoriented.nii.gz')
    mouse_annotation_path = os.path.join(MousePath, mouse_string + '_annotation.nii.gz')

    # Load images
    print('loading images')
    mouse_image = nib.load(mouse_path)
    mouse = mouse_image.get_fdata()
    mouse_annotation_image = nib.load(mouse_annotation_path)
    mouse_annotation = mouse_annotation_image.get_fdata()

    # Relevant variables to list
    mouse_image_list.append(mouse_image)
    mouse_list.append(mouse)
    mouse_annotation_image_list.append(mouse_annotation_image)
    mouse_annotation_list.append(mouse_annotation)

    # For each subject create perstructure folder if it not exists already
    perstructure_path = os.path.join(MousePath, 'perstructure')
    if not os.path.exists(perstructure_path):
            os.makedirs(perstructure_path)
            print(perstructure_path)

    # Preallocate 5D volumes with nan's
    mouse_defField_magnitude_5D = np.empty(list(reference_template_image.shape)+[3]+[nStructure])
    mouse_defField_magnitude_5D[:] = np.nan
    mouse_defField_magnitude_5D_list.append(mouse_defField_magnitude_5D)

# Loop through structures and isolate in reference and native space
# rigid flirt, affine flirt and syn isolated native structures to the isolated reference structure
# save relevant defField's
# superimpose inverse defField's with np.nansum(input_nanned)/np.sum(!np.isnan(input_nanned)) [averaging]
# save relevant superimosed defField' (reenter subject loop?)
reference_structure_crop_index_list = list()
reference_structure_mask_list = list()
for iVolume in range(reference_table.shape[0]):
    structure_name = reference_table.iloc[1]['name']
    structure_name = structure_name.replace(' ', '_')
    structure_name = structure_name.replace(',', '')

    qform = reference_template_image.affine
    children = reference_table.iloc[iVolume]['id_custom_children']
    reference_isolated_path = os.path.join(reference_path, 'perstructure', 'reference_template_test_' + structure_name + str(reference_table.iloc[iVolume]['id_custom']) + '.nii.gz')

    reference_isolated_image, reference_template_mask, reference_crop_index = \
        isolateCropStructure(reference_annotation, qform,
                             children,
                             reference_template,
                             reference_isolated_path)
    reference_isolated = reference_isolated_image.get_fdata()

    # Save cropping index so that later the inverted defField's to native space can be reassigned
    # for superposition with a 5D reference space volume filled with nan's
    # (filling with nan's done with mask)
    reference_structure_crop_index_list.append(reference_crop_index)
    reference_structure_mask_list.append(reference_template_mask)

    # Loop through mice
    for iMousePath, MousePath in enumerate(mouse_path_list):
        mouse_string = MousePath.split(os.sep)[-1]
        print(iMousePath)
        print(MousePath)
        print(datetime.datetime.now())

        # Define mouse paths
        mouse_image = mouse_image_list[iMousePath]
        mouse_array = mouse_list[iMousePath]
        mouse_annotation_image = mouse_annotation_image_list[iMousePath]
        mouse_annotation_array = mouse_annotation_list[iMousePath]
        native_isolated_path = os.path.join(MousePath, 'perstructure', mouse_string + '_isolated_' + structure_name + str(reference_table.iloc[iVolume]['id_custom']) + '.nii.gz')
        native_isolated_flirtRigid_path = os.path.join(MousePath, 'perstructure', mouse_string + '_isolated_flirtRigid_' + structure_name + str(reference_table.iloc[iVolume]['id_custom']) + '.mat')
        native_isolated_invflirtRigid_path = os.path.join(MousePath, 'perstructure', mouse_string + '_isolated_invflirtRigid_' + structure_name + str(reference_table.iloc[iVolume]['id_custom']) + '.mat')
        native_isolated_flirtedRigid_path = os.path.join(MousePath, 'perstructure', mouse_string + '_isolated_flirtedRigid_' + structure_name + str(reference_table.iloc[iVolume]['id_custom']) + '.nii.gz')
        native_isolated_flirt_path = os.path.join(MousePath, 'perstructure', mouse_string + '_isolated_flirt_' + structure_name + str(reference_table.iloc[iVolume]['id_custom']) + '.mat')
        native_isolated_invflirt_path = os.path.join(MousePath, 'perstructure', mouse_string + '_isolated_invflirt_' + structure_name + str(reference_table.iloc[iVolume]['id_custom']) + '.mat')
        native_isolated_flirted_path = os.path.join(MousePath, 'perstructure', mouse_string + '_isolated_flirted_' + structure_name + str(reference_table.iloc[iVolume]['id_custom']) + '.nii.gz')
        native_isolated_flirted_synned_path = os.path.join(MousePath, 'perstructure', mouse_string + '_isolated_flirted_synned_' + structure_name + str(reference_table.iloc[iVolume]['id_custom']) + '.nii.gz')

        # Crop and isolate
        mouse_isolated_image, mouse_isolated_mask, native_crop_index = \
            isolateCropStructure(mouse_annotation_array, mouse_image.affine,
                                 children,
                                 mouse_array,
                                 native_isolated_path)

        # FLIRT (rigid) cropped native isolate to cropped reference isolate
        print('FLIRT rigid start')
        os.system('flirt -in ' + native_isolated_path + ' \
                         -ref ' + reference_isolated_path + ' \
                         -out ' + native_isolated_flirtedRigid_path + ' \
                         -omat ' + native_isolated_flirtRigid_path + ' \
                         -dof ' + '6' + ' \
                         -verbose 0')    # FLIRT subject to reference
        os.system('convert_xfm -omat ' + native_isolated_invflirtRigid_path + ' -inverse ' + native_isolated_flirtRigid_path)

        # FLIRT (affine) cropped native isolate to cropped reference isolate
        print('FLIRT affine start')
        os.system('flirt -in ' + native_isolated_flirtedRigid_path + ' \
                         -ref ' + reference_isolated_path + ' \
                         -out ' + native_isolated_flirted_path + ' \
                         -omat ' + native_isolated_flirt_path + ' \
                         -verbose 0')
        os.system('convert_xfm -omat ' + native_isolated_invflirt_path + ' -inverse ' + native_isolated_flirt_path)
        native_isolated_flirted_image = nib.load(native_isolated_flirted_path)
        native_isolated_flirted = native_isolated_flirted_image.get_fdata()

        # SyN flirted images to reference
        print('SyN')
        level_iters = [10, 10, 10]
        CCMetric_radius = 4
        print(f'CCMetric_radius = {CCMetric_radius}')
        level_min_size = [round_half_up(num) for num in reference_isolated_image.shape/np.power(2, len(level_iters)-1)]
        while (2*CCMetric_radius+1) > np.min(np.array(level_min_size)):
            CCMetric_radius = CCMetric_radius - 1
            print(f'CCMetric_radius = {CCMetric_radius}')
        metric = CCMetric(3, radius=CCMetric_radius)
        sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)
        mapping = sdr.optimize(static=reference_isolated,
                               moving=native_isolated_flirted,
                               static_grid2world=reference_isolated_image.get_qform(),
                               moving_grid2world=native_isolated_flirted_image.get_qform())
        # with open(mouse_masked_flirted_syn_path, 'wb') as f:
        #     dump([mapping, metric, level_iters, sdr], f, protocol=4, compression='gzip')
        # with open(mouse_masked_syn_path, 'rb') as f:
        #     [mapping, metric, level_iters, sdr] = load(f)

        native_isolated_flirted_synned = mapping.transform(native_isolated_flirted)
        native_isolated_flirted_synned_image = nib.Nifti1Image(native_isolated_flirted_synned,
                                                               native_isolated_flirted_image.affine)
        nib.save(native_isolated_flirted_synned_image, native_isolated_flirted_synned_path)

        # Get FLIRT rigid vector field, assign to uncropped reference space, calculate magnitude and assign to 5D volume
        defField = imageFLIRT2defField(reference_template_image, native_isolated_invflirtRigid_path, dof='6')
        defField_magnitude = np.sqrt(np.sum(np.power(defField, 2), axis=3))
        mouse_defField_magnitude_5D_list[iMousePath][reference_crop_index[0]:defField_magnitude.shape[0],
                                                     reference_crop_index[1]:defField_magnitude.shape[1],
                                                     reference_crop_index[2]:defField_magnitude.shape[2],
                                                     0, iVolume] = defField_magnitude

        # Get FLIRT affine vector field, calculate magnitude and assign to 5D volume
        defField = imageFLIRT2defField(reference_template_image, native_isolated_invflirt_path)
        defField_magnitude = np.sqrt(np.sum(np.power(defField, 2), axis=3))
        mouse_defField_magnitude_5D_list[iMousePath][reference_crop_index[0]:defField_magnitude.shape[0],
                                                     reference_crop_index[1]:defField_magnitude.shape[1],
                                                     reference_crop_index[2]:defField_magnitude.shape[2],
                                                     1, iVolume] = defField_magnitude

        # Get SyN vector field, calculate magnitude and assign to 5D volume
        defField = mapping.get_backward_field()
        defField_magnitude = np.sqrt(np.sum(np.power(defField, 2), axis=3))
        mouse_defField_magnitude_5D_list[iMousePath][reference_crop_index[0]:defField_magnitude.shape[0],
                                                     reference_crop_index[1]:defField_magnitude.shape[1],
                                                     reference_crop_index[2]:defField_magnitude.shape[2],
                                                     2, iVolume] = defField_magnitude

# Loop again through subjects and average each of the three 4D volumes in the 5D volume list to get the final manual vector magnitudes
for iMousePath, MousePath in enumerate(mouse_path_list):
    mouse_string = MousePath.split(os.sep)[-1]
    native_isolated_flirtRigid_magnitude_path = os.path.join(MousePath, 'perstructure', mouse_string + '_flirtRigid_magnitude.nii.gz')
    native_isolated_flirt_magnitude_path = os.path.join(MousePath, 'perstructure', mouse_string + '_flirt_magnitude.nii.gz')
    native_isolated_syn_magnitude_path = os.path.join(MousePath, 'perstructure', mouse_string + '_syn_magnitude.nii.gz')

    defField_magnitude = np.squeeze(mouse_defField_magnitude_5D_list[iMousePath][:,:,:,0,:])
    nData = np.sum(np.logical_not(np.isnan(defField_magnitude)), axis=3)
    defField_magnitude_averaged = np.nansum(defField_magnitude)
    defField_magnitude_averaged = defField_magnitude_averaged/nData
    defField_magnitude_averaged[nData == 0] = 0
    save_image(defField_magnitude_averaged, reference_template_image, native_isolated_flirtRigid_magnitude_path)

    defField_magnitude = np.squeeze(mouse_defField_magnitude_5D_list[iMousePath][:,:,:,1,:])
    nData = np.sum(np.logical_not(np.isnan(defField_magnitude)), axis=3)
    defField_magnitude_averaged = np.nansum(defField_magnitude)
    defField_magnitude_averaged = defField_magnitude_averaged/nData
    defField_magnitude_averaged[nData == 0] = 0
    save_image(defField_magnitude_averaged, reference_template_image, native_isolated_flirt_magnitude_path)

    defField_magnitude = np.squeeze(mouse_defField_magnitude_5D_list[iMousePath][:,:,:,2,:])
    nData = np.sum(np.logical_not(np.isnan(defField_magnitude)), axis=3)
    defField_magnitude_averaged = np.nansum(defField_magnitude)
    defField_magnitude_averaged = defField_magnitude_averaged/nData
    defField_magnitude_averaged[nData == 0] = 0
    save_image(defField_magnitude_averaged, reference_template_image, native_isolated_syn_magnitude_path)
