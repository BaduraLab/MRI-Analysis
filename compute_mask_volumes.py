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
mouse_path = '/home/enzo/Desktop/Data/Mouse'
allen_path = '/usr/local/fsl/data/standard/allen_new'
data_path = '/home/enzo/Desktop/Data/Mouse/Processed_New'
analysis_path = '/home/enzo/Desktop/Data/Mouse/Analysis'
allen_image = os.path.join(allen_path, 'annotation_25_reoriented.nii.gz')
allen_image_flirted = os.path.join(allen_path, 'annotation_25_to_AMBMC_flirted.nii.gz')
picture_path = '/home/enzo/Desktop'
# voxel_volume = pow(0.05, 3)
voxel_volume = 0.000125
voxel_reference_volume = 1.5625e-5
voxel_AMBMC_volume = 3.375e-6

mouse_list = os.listdir(data_path)
nMouse = len(mouse_list)



## Mask volume computation, separate from Allen and invwarped computations
mouse_voxel_number = mouse_list.copy()
mouse_mask_volume = mouse_list.copy()
for iMouse, Mouse in enumerate(mouse_list):
    mouse_mask_image = nib.load(os.path.join(data_path, Mouse, (Mouse+'_mask_t=500_v=380_k=6.mask.nii.gz')))
    mouse_mask_image_array = mouse_mask_image.get_fdata()

    mouse_voxel_number[iMouse] = np.sum(mouse_mask_image_array>0)
    mouse_mask_volume[iMouse] = mouse_voxel_number[iMouse]*voxel_volume

mouse_mask_table = pd.DataFrame({'Mouse': mouse_list, 'MaskVoxelNumber': mouse_voxel_number, 'MaskVolume': mouse_mask_volume})
mouse_mask_table['Genotype'] = mouse_mask_table['Mouse'].str.split("_", n = 1, expand = True).iloc[:, 0]
mouse_mask_table['Sex'] = mouse_mask_table['Mouse'].str.split("_", n = 3, expand = True).iloc[:, 2]
mouse_mask_table.loc[mouse_mask_table['Sex'] != 'female', 'Sex'] = 'male'

# Plotting by genotype
gen_fig = plt.figure(1)
annotation_bin = nib.load(allen_image)
annotation_bin_array = annotation_bin.get_fdata()
allen_voxelnumber = np.sum(annotation_bin_array>.1)
allen_volume = allen_voxelnumber * voxel_reference_volume
mouse_mask_table_plot = pd.concat([mouse_mask_table, pd.DataFrame({'Mouse': ['Allen'], 'MaskVoxelNumber': [allen_voxelnumber], 'MaskVolume': [allen_volume], 'Genotype': ['Allen']})], sort=True)
# plt.suptitle('figure title')
# mouse_mask_table_plot=mouse_mask_table
ax = mouse_mask_table_plot[['MaskVolume', 'Genotype']].boxplot(by=['Genotype'])
# mouse_mask_table_plot.loc[mouse_mask_table_plot['Genotype'] == 'WT', 'Genotype'] = 1
# mouse_mask_table_plot.loc[mouse_mask_table_plot['Genotype'] == 'KO', 'Genotype'] = 2
plt.ylabel('$mm^3$')
plt.xlabel('Genotype')
plt.title('Brain mask volumes')
plt.suptitle('') # that's what you're after
# ax.set_xticklabels(['WT', 'KO'])
plt.show()
plt.savefig(os.path.join(analysis_path, 'Boxplot_MaskVolumes_ByGenotype_Allen'))

# Plotting by genotype and sex
gensex_fig = plt.figure(2)
mouse_mask_table[['MaskVolume', 'Genotype', 'Sex']].boxplot(by=['Genotype', 'Sex'])
plt.ylabel('$mm^3$')
plt.ylabel('')
plt.savefig(os.path.join(analysis_path, 'Boxplot_MaskVolumes_ByGenotypeSex'))
# plt.show()

# pval calculation, equal_var for now ;)
mouse_mask_table[mouse_mask_table['Genotype'] == 'WT']['MaskVolume']
cat2 = mouse_mask_table[mouse_mask_table['Genotype'] == 'KO']

print(ttest_ind(mouse_mask_table[mouse_mask_table['Genotype'] == 'WT']['MaskVolume'],
                mouse_mask_table[mouse_mask_table['Genotype'] == 'KO']['MaskVolume'],
                equal_var=True))

mouse_mask_table.to_csv(os.path.join(analysis_path, 'Mouse_maskvolume_table.csv'))





## Function to compute volumes for image
def image2volumetable(image_path, voxel_volume):
    # Compute voxel numbers and volumes and output to table
    mouse_mask_image = nib.load(image_path)
    mouse_mask_image_array = mouse_mask_image.get_fdata()
    [mouse_volume_integer, mouse_voxel_number] = np.unique(np.int64(np.round(mouse_mask_image_array)),
                                                           return_counts=True)
    mouse_volume = mouse_voxel_number * voxel_volume
    mouse_table_reference = pd.DataFrame(
        {'Mouse': 'allen', 'VolumeInteger': mouse_volume_integer, 'VoxelNumber': mouse_voxel_number,
         'Volume': mouse_volume})

    # print number of nonzero volumes
    print(mouse_table_reference.shape[0])
    print(np.sum(mouse_table_reference['Volume'] != 0))
    print(np.sum(mouse_table_reference['Volume'][mouse_table_reference['VolumeInteger']!=0]))

    # Attach parent path to volumes to table - volumes which are not in structure graph are removed
    mouse_table_reference = pd.merge(left=structure_graph, right=mouse_table_reference,
                                     left_on='id_custom', right_on='VolumeInteger', how='outer')
    mouse_table_reference.loc[np.isnan(mouse_table_reference['VoxelNumber']), 'VoxelNumber'] = 0
    mouse_table_reference.loc[np.isnan(mouse_table_reference['Volume']), 'Volume'] = 0

    # print number of nonzero volumes remaining after merge, which should be the same
    print(np.sum(mouse_table_reference['Volume'] != 0))

    # Fill in structures without a volume by summing relevant lower level structures
    for iVolume in range(mouse_table_reference.shape[0]):
        # only include structures with volumes here
        if (not np.isnan(mouse_table_reference.loc[iVolume, 'VolumeInteger'])) & isinstance(
                mouse_table_reference.loc[iVolume, 'structure_id_path_custom'], str):
            for iParent in list(map(int, mouse_table_reference.loc[iVolume, 'structure_id_path_custom'].strip('][').split(', ')))[:-1]:
                # Add volume to each parent
                mouse_table_reference.loc[mouse_table_reference['id_custom'] == iParent, 'Volume'] += \
                mouse_table_reference.loc[iVolume, 'Volume']
                mouse_table_reference.loc[mouse_table_reference['id_custom'] == iParent, 'VoxelNumber'] += \
                mouse_table_reference.loc[iVolume, 'VoxelNumber']

    return mouse_table_reference

# Define structure graph
structure_graph = pd.read_csv(os.path.join(allen_path, 'structure_graph_remapped.csv'))



## Get list of invwarped annotation files and compute volumes for them, compiling eventually into single table
mouse_invwarped_list = glob.glob(data_path+'/*/*/*invwarped*')
mouse_table_invwarped_list = mouse_invwarped_list.copy()
for iMouse, Mouse in enumerate(mouse_invwarped_list):
    Mouse_name = Mouse.split('/')[7]
    FNIRT_run = Mouse.split('/')[-2]

    mouse_table_invwarped_list[iMouse] = image2volumetable(Mouse, voxel_volume)
    mouse_table_invwarped_list[iMouse]['Mouse'] = Mouse_name
    mouse_table_invwarped_list[iMouse]['FNIRT_run'] = FNIRT_run

mouse_table_invwarped = pd.concat(mouse_table_invwarped_list)

mouse_mask_table.to_csv(os.path.join(analysis_path, 'Mouse_invwarped_table.csv'))



## Compute volumes for allen reference
allen_table = image2volumetable(allen_image, voxel_reference_volume)
allen_table['Mouse'] = 'allen'
allen_table.to_csv(os.path.join(analysis_path, 'allen_table.csv'))



## Compute volumes for flirted allen reference
allen_flirted_table = image2volumetable(allen_image_flirted, voxel_AMBMC_volume)
allen_flirted_table['Mouse'] = 'allen_flirted'
allen_flirted_table.to_csv(os.path.join(analysis_path, 'allen_flirted_table.csv'))






# for i in range(100):
#     print(i)
#     print(np.any(pd.DataFrame(oapi.get_structures(i+1))['id']==182305696))
