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
data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
reference_path = os.path.join('Data', 'Mouse', 'Reference')
analysis_path = os.path.join('Data', 'Mouse', 'Analysis')
mouse_path_list = glob.glob(os.path.join(data_path, '*', '*invsynned*cerebellum.nii.gz'))
reference_structure_path = os.path.join(reference_path, 'structure_graph_remapped_lowdetail.csv')
analysis_table_path = os.path.join(analysis_path, 'all_volumes.csv')
# VOIs = ['Lobule II', 'Lobules IV-V', 'Substantia nigra, compact part', 'Substantia nigra, reticular part']
# VOIs = ['Substantia nigra, compact part', 'Substantia nigra, reticular part']


# Load
structure = pd.read_csv(reference_structure_path)
volume_table = pd.read_csv(analysis_table_path)
VOIs = volume_table['name']
VOIs = np.unique(VOIs[np.logical_not(pd.isnull(VOIs))])
# volume_table = pd.merge(left=volume_table, right=structure.loc[:, ['name', 'id_custom']],
#                          left_on='name', right_on='name')







# Go through mice and compute center of mass of VOIs
mouse_path_list = [os.path.join(reference_path, 'annotation_50_reoriented.nii.gz')]
structure_centers_table_list = list()
for iMousePath, MousePath in enumerate(mouse_path_list):
    mouse_string = 'allen'
    # mouse_string = MousePath.split(os.sep)[-2]
    volume_table_mouse = volume_table[volume_table['Mouse'] == 'WT_50']

    # Load image
    mouse_image = nib.load(MousePath)
    mouse = mouse_image.get_fdata()

    print(mouse_string)

    # Calculate center of mouse brain mask, use y value to separate hemispheres
    mask_center = ndimage.measurements.center_of_mass(mouse>0)

    # Loop over VOIs
    mouse_arr = [None for i in range(len(VOIs))]
    name_arr = mouse_arr.copy()
    CoM_right_arr = mouse_arr.copy()
    CoM_left_arr = mouse_arr.copy()
    CoM_arr = mouse_arr.copy()
    for iVOI, VOI in enumerate(VOIs):
        # Get id_custom of VOI
        mouse_VOI_path = os.path.join(analysis_path, mouse_string + '_' + VOI)
        volume_table_mouse_structure = volume_table_mouse[volume_table_mouse['name'] == VOI]
        VOI_id_custom = volume_table_mouse_structure['id_custom'].iloc[0]

        # for each VOI, determine its center of mass in both affine and voxel coordinates
        mouse_VOI_logical = mouse == VOI_id_custom
        VOI_indices = np.array(np.where(mouse_VOI_logical))

        right_logical = VOI_indices[0]>mask_center[0]
        left_logical = np.logical_not(right_logical)

        VOI_indices_center_right = np.mean(VOI_indices[:, right_logical], 1)
        VOI_indices_center_left = np.mean(VOI_indices[:, left_logical], 1)
        VOI_indices_center = ndimage.measurements.center_of_mass(mouse_VOI_logical)
        VOI_indices_center_affine_right = mouse_image.affine.dot(np.concatenate((VOI_indices_center_right, np.array([0]))))
        VOI_indices_center_affine_left = mouse_image.affine.dot(np.concatenate((VOI_indices_center_left, np.array([0]))))
        VOI_indices_center_affine = mouse_image.affine.dot(np.concatenate((VOI_indices_center, np.array([0]))))

        # print(mouse_string)
        # print(VOI)
        # print(VOI_indices_center_right)
        # print(VOI_indices_center_left)
        # print(VOI_indices_center)
        # # print(VOI_indices_center2)
        # print(VOI_indices_center_affine_right)
        # print(VOI_indices_center_affine_left)
        # print(VOI_indices_center_affine)

        mouse_arr[iVOI] = mouse_string
        name_arr[iVOI] = VOI
        CoM_right_arr[iVOI] = VOI_indices_center_right
        CoM_left_arr[iVOI] = VOI_indices_center_left
        CoM_arr[iVOI] = VOI_indices_center

    mouse_structure_table = pd.DataFrame({'mouse': mouse_arr,
                                          'name': name_arr,
                                          'CoM': CoM_arr,
                                          'CoM_right': CoM_right_arr,
                                          'CoM_left': CoM_left_arr})
    structure_centers_table_list.append(mouse_structure_table)


structure_centers_table = pd.concat(structure_centers_table_list)
structure_centers_table.to_csv(os.path.join(analysis_path, 'reference_CoM.csv'))



# reference_CoM = pd.read_csv(os.path.join(analysis_path, 'reference_CoM.csv'))
reference_CoM = structure_centers_table
reference_CoM_filt = reference_CoM[~(np.array([np.any(np.isnan(reference_CoM['CoM_right'][i])) for i in range(len(reference_CoM['CoM_right']))]) |
                                     np.array([np.any(np.isnan(reference_CoM['CoM_left'][i])) for i in range(len(reference_CoM['CoM_left']))]) )]
reference_CoM_D = (dist.pairwise(np.vstack(np.array(reference_CoM_filt['CoM_right'])))+dist.pairwise(np.vstack(np.array(reference_CoM_filt['CoM_left']))))/2

reference_CoM_invD = np.log(reference_CoM_D)
np.fill_diagonal(reference_CoM_invD, 0)
reference_CoM_invD = reference_CoM_invD - np.min(reference_CoM_invD)
reference_CoM_invD = 1 - (reference_CoM_invD / np.max(reference_CoM_invD))
np.fill_diagonal(reference_CoM_invD, 0)

structs_1D_space = embedding.fit_transform(reference_CoM_invD)
structs_1D_space_int = np.argsort(np.squeeze(structs_1D_space))
structs_1D_space_int = np.round(structs_1D_space_int/np.max(structs_1D_space_int)*(len(VOIs)-200-1)).astype('int')+1+200
structs_1D_space_int_rest = np.setxor1d(np.arange(len(VOIs))+1, np.sort(structs_1D_space_int))
np.random.shuffle(structs_1D_space_int_rest)
# everything that does not have an integer should be assigned an integer according to rest
# after every structure has an integer, do a remap with special numpy
# save this remapping also to the csv ("_MDS"?)
reference_CoM_filt.insert(2, 'id_MDS', structs_1D_space_int)
structure_MDS = pd.merge(left=structure, left_on='name',
         right=reference_CoM_filt[['name', 'id_MDS']], right_on='name',
         how='left')
structure_MDS.loc[pd.isnull(structure_MDS['id_MDS']), 'id_MDS'] = structs_1D_space_int_rest
structure_MDS = structure_MDS.astype({'id_MDS': int})
structure_MDS.to_csv(os.path.join(reference_path, 'structure_graph_MDS.csv'))

mouse_path_list = glob.glob(os.path.join(data_path, '*', '*invsynned*cerebellum.nii.gz')) + \
    [os.path.join(reference_path, 'annotation_50_reoriented.nii.gz')]
for Mouse in mouse_path_list:
    print(Mouse)

    output_path = Mouse.split('.')[0]+'_MDS.nii.gz'

    input_image = nib.load(Mouse)
    input = input_image.get_fdata()

    input_shape = input.shape
    input = input.reshape(-1)
    input = npi.remap(input, list(structure_MDS['id_custom']), list(structure_MDS['id_MDS']))
    input = input.reshape(input_shape)
    output_image = nib.Nifti1Image(input, input_image.affine, input_image.header)
    nib.save(output_image, output_path)





# for Mouse in np.unique(structure_centers_table['mouse']):
    # structure_centers_mouse_table = structure_centers_table[structure_centers_table['mouse']==Mouse]



        # # for each VOI, isolate it and create new nifti files
        # mouse_VOI = mouse.copy()
        # mouse_VOI[np.logical_not(mouse_VOI_logical)] = 0
        # mouse_VOI_image = nib.Nifti1Image(mouse_VOI, mouse_image.affine, mouse_image.header)
        # nib.save(mouse_VOI_image, mouse_VOI_path)
        # # view in mricrogl outside of python


# volume_table[np.logical_and(np.isin(volume_table['name'], VOIs), np.isin(volume_table['Mouse'], ['WT_50', 'KO_6']))]

# # Plot boxplot plots for VOIs
# volume_name = 'Lobule II'
# ax = mouse_table_all[mouse_table_all['name']==volume_name][['Volume', 'Genotype']].boxplot(by=['Genotype'])
# plt.ylabel('$mm^3$')
# plt.xlabel('Genotype')
# plt.title(volume_name + ' volumes')
# plt.suptitle('') # that's what you're after
# # ax.set_xticklabels(['WT', 'KO'])
# # plt.show()
# plt.savefig(os.path.join(analysis_path, 'Boxplot_'+volume_name+'_ByGenotype'))
#
# volume_name = 'Substantia nigra, compact part'
# ax = mouse_table_all[mouse_table_all['name']==volume_name][['Volume', 'Genotype']].boxplot(by=['Genotype'])
# plt.ylabel('$mm^3$')
# plt.xlabel('Genotype')
# plt.title(volume_name + ' volumes')
# plt.suptitle('') # that's what you're after
# # ax.set_xticklabels(['WT', 'KO'])
# # plt.show()
# plt.savefig(os.path.join(analysis_path, 'Boxplot_'+volume_name+'_ByGenotype'))
#
# volume_name = 'Substantia nigra, reticular part'
# ax = mouse_table_all[mouse_table_all['name']==volume_name][['Volume', 'Genotype']].boxplot(by=['Genotype'])
# plt.ylabel('$mm^3$')
# plt.xlabel('Genotype')
# plt.title(volume_name + ' volumes')
# plt.suptitle('') # that's what you're after
# # ax.set_xticklabels(['WT', 'KO'])
# # plt.show()
# plt.savefig(os.path.join(analysis_path, 'Boxplot_'+volume_name+'_ByGenotype'))
#
# # Plotting by genotype and sex
# volume_name = 'Lobules IV-V'
# ax = mouse_table_all[mouse_table_all['name']==volume_name][['Volume', 'Genotype', 'Sex']].boxplot(by=['Genotype', 'Sex'])
# plt.ylabel('$mm^3$')
# plt.xlabel('Genotype and Sex')
# plt.title(volume_name + ' volumes')
# plt.suptitle('') # that's what you're after
# # ax.set_xticklabels(['WT', 'KO'])
# # plt.show()
# plt.savefig(os.path.join(analysis_path, 'Boxplot_'+volume_name+'_ByGenotypeSex'))