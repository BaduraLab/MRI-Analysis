import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt
import glob
import csv
from scipy.stats import ttest_ind



# Define
mouse_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse'
allen_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse/allen'
data_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse/Processed_New'
analysis_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse/Analysis'
allen_image = os.path.join(allen_path, 'annotation_25_reoriented.nii.gz')
allen_image_flirted = os.path.join(allen_path, 'annotation_25_to_AMBMC_flirted.nii.gz')
# voxel_volume = pow(0.05, 3)
voxel_volume = 0.000125
voxel_reference_volume = 1.5625e-5

mouse_list = os.listdir(data_path)
nMouse = len(mouse_list)



mouse_voxel_number = mouse_list.copy()
mouse_mask_volume = mouse_list.copy()
for iMouse, Mouse in enumerate(mouse_list):
    mouse_mask_image = nib.load(os.path.join(data_path, Mouse, (Mouse+'_mask_t=500_v=380_k=6.mask.nii.gz')))
    mouse_mask_image_array = mouse_mask_image.get_fdata()

    mouse_voxel_number[iMouse] = np.sum(mouse_mask_image_array>0)
    mouse_mask_volume[iMouse] = mouse_voxel_number[iMouse]*voxel_volume

mouse_table = pd.DataFrame({'Mouse': mouse_list, 'MaskVoxelNumber': mouse_voxel_number, 'MaskVolume': mouse_mask_volume})
mouse_table['Genotype'] = mouse_table['Mouse'].str.split("_", n = 1, expand = True).iloc[:, 0]
mouse_table['Sex'] = mouse_table['Mouse'].str.split("_", n = 3, expand = True).iloc[:, 2]
mouse_table.loc[mouse_table['Sex']!='female', 'Sex'] = 'male'



plt.figure()
mouse_table[['MaskVolume', 'Genotype']].boxplot(by=['Genotype'])
plt.savefig('C:/Users/enzo/Downloads/Study/Current Courses/MEP/Pictures/Boxplot_MaskVolumes_ByGenotype')
plt.ylabel('$mm^3$')
plt.show()

plt.figure()
mouse_table[['MaskVolume', 'Genotype', 'Sex']].boxplot(by=['Genotype', 'Sex'])
plt.savefig('C:/Users/enzo/Downloads/Study/Current Courses/MEP/Pictures/Boxplot_MaskVolumes_ByGenotypeSex')
plt.ylabel('$mm^3$')
plt.show()



mouse_table[mouse_table['Genotype'] == 'WT']['MaskVolume']
cat2 = mouse_table[mouse_table['Genotype'] == 'KO']

print(ttest_ind(mouse_table[mouse_table['Genotype'] == 'WT']['MaskVolume'],
                mouse_table[mouse_table['Genotype'] == 'KO']['MaskVolume'],
                equal_var=True))

## Get list of invwarped annotation files
example_invwarped = os.path.join(data_path, '')
mouse_invwarped_list = glob.glob(data_path+'/*/*/*invwarped*')

mouse_volume_integer = mouse_invwarped_list.copy()
mouse_voxel_number = mouse_invwarped_list.copy()
mouse_mask_volume = mouse_invwarped_list.copy()
mouse_table_invwarped = mouse_invwarped_list.copy()
for iMouse, Mouse in enumerate(mouse_invwarped_list):
    Mouse_name = Mouse.split('\\')[1]

    mouse_mask_image = nib.load(Mouse)
    mouse_mask_image_array = mouse_mask_image.get_fdata()

    [mouse_volume_integer[iMouse], mouse_voxel_number[iMouse]] = np.unique(np.int64(np.round(mouse_mask_image_array)), return_counts=True)
    mouse_mask_volume[iMouse] = mouse_voxel_number[iMouse] * voxel_volume

    mouse_table_invwarped[iMouse] = pd.DataFrame({'Mouse': Mouse_name, 'VolumeInteger': mouse_volume_integer[iMouse], 'VoxelNumber': mouse_voxel_number[iMouse], 'Volume': mouse_mask_volume[iMouse]})

mouse_table_invwarped = pd.concat(mouse_table_invwarped)

structure_graph = pd.read_csv('C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse/allen/structure_graph.csv')

mouse_table_invwarped_2 = pd.merge(left=structure_graph, right=mouse_table_invwarped,
                                   left_on='id', right_on='VolumeInteger', how='left')
mouse_table_invwarped_2.loc[np.isnan(mouse_table_invwarped_2['VoxelNumber']), 'VoxelNumber'] = 0
mouse_table_invwarped_2.loc[np.isnan(mouse_table_invwarped_2['Volume']), 'Volume'] = 0

## Compute Allen Atlas volumes for comparison
mouse_mask_image = nib.load(allen_image)
mouse_mask_image_array = mouse_mask_image.get_fdata()

[mouse_volume_integer, mouse_voxel_number] = np.unique(np.int64(np.round(mouse_mask_image_array)),
                                                       return_counts=True)
mouse_mask_volume = mouse_voxel_number * voxel_reference_volume

mouse_table_reference = pd.DataFrame(
    {'Mouse': 'allen', 'VolumeInteger': mouse_volume_integer, 'VoxelNumber': mouse_voxel_number,
     'Volume': mouse_mask_volume})



## Compute Allen Atlas volumes for comparison

def image2volumetable(image_path, voxel_volume):
    mouse_mask_image = nib.load(allen_image_flirted)
    mouse_mask_image_array = mouse_mask_image.get_fdata()

    [mouse_volume_integer, mouse_voxel_number] = np.unique(np.int64(np.round(mouse_mask_image_array)),
                                                           return_counts=True)
    mouse_mask_volume = mouse_voxel_number * voxel_reference_volume

    mouse_table_reference = pd.DataFrame(
        {'Mouse': 'allen', 'VolumeInteger': mouse_volume_integer, 'VoxelNumber': mouse_voxel_number,
         'Volume': mouse_mask_volume})

    # loop over all structures
    for iVolume in range(mouse_table_reference.shape[0]):
        # only include structures with volumes here
        if not np.isnan(mouse_table_reference.loc[iVolume, 'VolumeInteger']):
            for iParent in list(map(int, mouse_table_reference.loc[iVolume, 'structure_id_path'].strip('][').split(', ')))[:-1]:
                # Add volume to each parent
                mouse_table_reference.loc[iParent, 'Volume'] += mouse_table_reference.loc[iVolume, 'Volume']
                mouse_table_reference.loc[iParent, 'VoxelNumber'] += mouse_table_reference.loc[iVolume, 'VoxelNumber']



## Write tables
mouse_table.to_csv(os.path.join(analysis_path, 'Mouse_maskvolume_table.csv'))
mouse_table_invwarped_2.to_csv(os.path.join(analysis_path, 'Mouse_volume_table.csv'))

# add mask volumes to warped volumes?
# add reference volumes for comparison, not voxel_volume difference
# function to make volume table for mouse/reference input voxel_volume and image_path, output table with filled in volumes