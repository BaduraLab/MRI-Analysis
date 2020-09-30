import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import glob
import csv
from scipy.stats import ttest_ind
from scipy import ndimage
from sklearn.neighbors import DistanceMetric
dist = DistanceMetric.get_metric('euclidean')
from sklearn.manifold import MDS
embedding = MDS(n_components=1, dissimilarity='precomputed')
from pathlib import Path
import numpy_indexed as npi



# Define
data_path = os.path.join('Data', 'Human', 'Processed')
reference_path = os.path.join('Data', 'Human', 'Reference')
analysis_path = os.path.join('Data', 'Human', 'Analysis')
structure_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT.txt'),
                       os.path.join(reference_path, 'subcortical', 'subcortical.csv'),
                       os.path.join(reference_path, 'CerebrA', 'CerebrA.csv')]
structure_MDS_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT_MDS.csv'),
                       os.path.join(reference_path, 'subcortical', 'subcortical_MDS.csv'),
                       os.path.join(reference_path, 'CerebrA', 'CerebrA_MDS.csv')]
# analysis_table_path = os.path.join(analysis_path, 'all_volumes.csv')
structure_table_list = [pd.read_csv(structure_path_list[0],
                                    delim_whitespace=True,
                                    names=['VolumeInteger', 'name', 'ID']),
                        pd.read_csv(structure_path_list[1],
                                    names=['VolumeInteger', 'acronym', 'name']),
                        pd.read_csv(structure_path_list[2])]
# VOIs = ['Lobule II', 'Lobules IV-V', 'Substantia nigra, compact part', 'Substantia nigra, reticular part']
# VOIs = ['Substantia nigra, compact part', 'Substantia nigra, reticular part']
annotation_path_list = [os.path.join(reference_path,
                                     'suit',
                                     'atlasesSUIT',
                                     'Lobules-SUIT.nii'),
                        os.path.join(reference_path,
                                     'subcortical',
                                     'prob_atlas_bilateral_thrarg_0.4.nii.gz'),
                        os.path.join(reference_path,
                                     'CerebrA',
                                     'mni_icbm152_CerebrA_tal_nlin_sym_09c_reoriented.nii.gz')]
annotation_name_list = ['suit',
                        'subcortical',
                        'CerebrA']





# Go through annotation files (reference thrarg) and compute center of mass of VOIs
CoM_structure_table_list = list()
for iAnnotationPath, AnnotationPath in enumerate(annotation_path_list):
    annotation_string = annotation_name_list[iAnnotationPath]
    # mouse_string = MousePath.split(os.sep)[-2]
    structure_table = structure_table_list[iAnnotationPath]
    VOIs = structure_table['name']

    # Load image
    annotation_image = nib.load(AnnotationPath)
    annotation = annotation_image.get_fdata()

    print(annotation_string)

    # Calculate center of mouse brain mask, use y value to separate hemispheres
    mask_center = ndimage.measurements.center_of_mass(annotation > 0)

    # Loop over VOIs
    annotation_arr = [None for i in range(len(VOIs))]
    name_arr = annotation_arr.copy()
    VolumeInteger_arr = annotation_arr.copy()
    CoM_right_arr = annotation_arr.copy()
    CoM_left_arr = annotation_arr.copy()
    CoM_arr = annotation_arr.copy()
    for iVOI, VOI in enumerate(VOIs):
        # Get id_custom of VOI
        annotation_VOI_path = os.path.join(analysis_path, annotation_string + '_' + VOI)
        structure_VOI_table = structure_table[VOIs == VOI]
        VOI_id_custom = structure_VOI_table['VolumeInteger'].iloc[0]

        # for each VOI, determine its center of mass in both affine and voxel coordinates
        mouse_VOI_logical = annotation == VOI_id_custom
        VOI_indices = np.array(np.where(mouse_VOI_logical))

        if iAnnotationPath == 1:
            right_logical = VOI_indices[0]>mask_center[0]
            left_logical = np.logical_not(right_logical)

            VOI_indices_center_right = np.mean(VOI_indices[:, right_logical], 1)
            VOI_indices_center_left = np.mean(VOI_indices[:, left_logical], 1)
            VOI_indices_center = ndimage.measurements.center_of_mass(mouse_VOI_logical)
            VOI_indices_center_affine_right = annotation_image.affine.dot(np.concatenate((VOI_indices_center_right, np.array([0]))))
            VOI_indices_center_affine_left = annotation_image.affine.dot(np.concatenate((VOI_indices_center_left, np.array([0]))))
            VOI_indices_center_affine = annotation_image.affine.dot(np.concatenate((VOI_indices_center, np.array([0]))))

            annotation_arr[iVOI] = annotation_string
            name_arr[iVOI] = VOI
            VolumeInteger_arr[iVOI] = VOI_id_custom
            CoM_right_arr[iVOI] = VOI_indices_center_right
            CoM_left_arr[iVOI] = VOI_indices_center_left
            CoM_arr[iVOI] = VOI_indices_center
        else:
            VOI_indices_center = ndimage.measurements.center_of_mass(mouse_VOI_logical)
            VOI_indices_center_affine = annotation_image.affine.dot(np.concatenate((VOI_indices_center, np.array([0]))))

            annotation_arr[iVOI] = annotation_string
            name_arr[iVOI] = VOI
            VolumeInteger_arr[iVOI] = VOI_id_custom
            CoM_arr[iVOI] = VOI_indices_center


        # print(mouse_string)
        # print(VOI)
        # print(VOI_indices_center_right)
        # print(VOI_indices_center_left)
        # print(VOI_indices_center)
        # # print(VOI_indices_center2)
        # print(VOI_indices_center_affine_right)
        # print(VOI_indices_center_affine_left)
        # print(VOI_indices_center_affine)

    if iAnnotationPath == 1:
        CoM_structure_table = pd.DataFrame({'annotation': annotation_arr,
                                            'name': name_arr,
                                            'VolumeInteger': VolumeInteger_arr,
                                            'CoM': CoM_arr,
                                            'CoM_right': CoM_right_arr,
                                            'CoM_left': CoM_left_arr})
        CoM_structure_D = (dist.pairwise(np.vstack(np.array(CoM_structure_table['CoM_right']))) +
                           dist.pairwise(np.vstack(np.array(CoM_structure_table['CoM_left'])))) / 2
    else:
        CoM_structure_table = pd.DataFrame({'annotation': annotation_arr,
                                            'name': name_arr,
                                            'VolumeInteger': VolumeInteger_arr,
                                            'CoM': CoM_arr})
        CoM_structure_D = dist.pairwise(np.vstack(np.array(CoM_structure_table['CoM'])))

    CoM_structure_invD = np.log(CoM_structure_D)
    np.fill_diagonal(CoM_structure_invD, 0)
    CoM_structure_invD = CoM_structure_invD - np.min(CoM_structure_invD)
    CoM_structure_invD = 1 - (CoM_structure_invD / np.max(CoM_structure_invD))
    np.fill_diagonal(CoM_structure_invD, 0)

    structs_1D_space = embedding.fit_transform(CoM_structure_invD)
    structs_1D_space_int = np.argsort(np.squeeze(structs_1D_space)) + 1
    CoM_structure_table.insert(3, 'VolumeInteger_MDS', structs_1D_space_int)
    CoM_structure_table.astype({'VolumeInteger_MDS': int})
    CoM_structure_table.to_csv(structure_MDS_path_list[iAnnotationPath])

    CoM_structure_table_list.append(CoM_structure_table)

CoM_structure_table_all = pd.concat(CoM_structure_table_list, sort=True)
CoM_structure_table_all.to_csv(os.path.join(analysis_path, 'reference_CoM.csv'))



# Apply MDS on human output (use CoM_structure_table_list)
input_path_list_list = [glob.glob(os.path.join(data_path, '*', '*annotation_orsuit_thrarg_adjusted.nii.gz')) + [annotation_path_list[0]],
                        glob.glob(os.path.join(data_path, '*', '*annotation_subcortical_thrarg.nii.gz')) + [annotation_path_list[1]],
                        glob.glob(os.path.join(data_path, '*', '*annotation_CerebrA_thrarg_adjusted.nii.gz')) + [annotation_path_list[2]]]
for iInput, input_path_list in enumerate(input_path_list_list):
    print(iInput)
    CoM_structure_table = CoM_structure_table_list[iInput]

    for InputPath in input_path_list:
        output_path = InputPath.split('.')[0]+'_MDS.nii.gz'
        print(output_path)

        input_image = nib.load(InputPath)
        input = input_image.get_fdata()

        input = np.round(input).astype(int) # ensure integer annotation input

        input_shape = input.shape
        input = input.reshape(-1)
        print(list(CoM_structure_table['VolumeInteger']))
        print(list(CoM_structure_table['VolumeInteger_MDS']))
        input = npi.remap(input, list(CoM_structure_table['VolumeInteger']), list(CoM_structure_table['VolumeInteger_MDS']))
        input = input.reshape(input_shape)
        output_image = nib.Nifti1Image(input, input_image.affine, input_image.header)
        nib.save(output_image, output_path)
