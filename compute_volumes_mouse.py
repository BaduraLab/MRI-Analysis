import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import glob
import csv
from scipy.stats import ttest_ind
from pathlib import Path
import numpy_indexed as npi
from functions import remap_3D
from functions import save_image



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
    mouse_table_reference = pd.merge(left=structure, right=mouse_table_reference,
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



# Define
data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
reference_path = os.path.join('Data', 'Mouse', 'Reference')
analysis_path = os.path.join('Data', 'Mouse', 'Analysis')
mouse_path_list = glob.glob(os.path.join(data_path, '*', '*invsynned*cerebellum_lobular.nii.gz'))
reference_structure_path = os.path.join(reference_path, 'structure_graph_remapped_lowdetail.csv')
voxel_reference_volume = 0.000125
voxel_volume = 0.000125
annotation_path = os.path.join(reference_path, 'annotation_50_reoriented.nii.gz')

# Follows
nMouse = len(mouse_path_list)
structure = pd.read_csv(reference_structure_path)
Path(os.path.join(analysis_path, 'perstructure')).mkdir(parents=True, exist_ok=True)



# Reference volumes
reference_table = image2volumetable(annotation_path, voxel_reference_volume)
reference_table.to_csv(os.path.join(analysis_path, 'reference_volumes.csv'))

reference_table['in_cerebellum'] = False
for iVolume in range(reference_table.shape[0]):
    if isinstance(reference_table.loc[iVolume, 'structure_id_path'], str):
        reference_table.loc[iVolume, 'in_cerebellum'] = 512 in list(map(int, reference_table.loc[iVolume, 'structure_id_path'].strip('][').split(', ')))
reference_cerebellum_table = reference_table[reference_table['in_cerebellum']][['name', 'acronym', 'id', 'structure_id_path', 'id_custom', 'structure_id_path_custom', 'VoxelNumber', 'Volume']]
reference_cerebellum_table.to_csv(os.path.join(analysis_path, 'reference_volumes_cerebellum.csv'))



mouse_table_list = list()
for iMouse, Mouse in enumerate(mouse_path_list):
    subject = Mouse.split(os.path.sep)[-2]
    print(Mouse)
    print(subject)
    mouse_table = image2volumetable(Mouse, voxel_volume)
    mouse_table.to_csv(os.path.join(analysis_path, subject+'_volumes.csv'))

    mouse_table['Mouse'] = subject
    mouse_table['Genotype'] = subject.split('_')[0]
    mouse_table['Sex'] = subject.split('_')[-1]
    mouse_table.loc[mouse_table['Sex']!='female', 'Sex'] = 'male'
    mouse_table = mouse_table[['Mouse', 'Genotype', 'Sex', 'name', 'acronym', 'id_custom', 'structure_id_path_custom',
                               'VoxelNumber', 'Volume']]

    mouse_table_list.append(mouse_table)

mouse_table_all = pd.concat(mouse_table_list, ignore_index=True)

reference_table['rVolume'] = reference_table['Volume']
output_table_all = pd.merge(left=mouse_table_all,
                            right=reference_table.loc[:, ['id_custom', 'name', 'rVolume']],
                            left_on=['id_custom', 'name'],
                            right_on=['id_custom', 'name'])
output_table_all['VolumeNormalized'] = output_table_all['Volume']/output_table_all['rVolume']

output_table_all.to_csv(os.path.join(analysis_path, 'all_volumes.csv'))

output_table_cerebellum = pd.merge(left=output_table_all, right=reference_table[['id_custom', 'in_cerebellum']],
         left_on='id_custom', right_on='id_custom')
output_table_cerebellum = output_table_cerebellum[output_table_cerebellum['in_cerebellum']]
output_table_cerebellum = output_table_cerebellum.drop('in_cerebellum', 1)
output_table_cerebellum.to_csv(os.path.join(analysis_path, 'cerebellum_volumes.csv'))


# In volume table, go through each structure and determine the p-value between genotypes, create a new p-value table
mouse_table_pername_list = list()
mouse_table_all_nobackground = mouse_table_all.loc[np.logical_not(pd.isnull(mouse_table_all['name']))]
# mouse_table_all.loc[pd.isnull(mouse_table_all['name']), 'name'] = 'background'
for nameStruct in np.unique(np.array(mouse_table_all_nobackground['name'].astype('category'))):
    mouse_table_nameStruct = mouse_table_all_nobackground[mouse_table_all_nobackground['name'] == nameStruct]
    mouse_table_nameStruct_WT = mouse_table_nameStruct.loc[mouse_table_all_nobackground['Genotype'] == 'WT']
    mouse_table_nameStruct_KO = mouse_table_nameStruct.loc[mouse_table_all_nobackground['Genotype'] == 'KO']
    [t_stat, p_val] = ttest_ind(mouse_table_nameStruct_WT['Volume'],
        mouse_table_nameStruct_KO['Volume'],
        equal_var=False)
    mean_WT = np.mean(mouse_table_nameStruct_WT['Volume'])
    mean_KO = np.mean(mouse_table_nameStruct_KO['Volume'])
    mouse_table_pername_list.append(pd.DataFrame({'name': [nameStruct],
                                                  'pVal': [p_val],
                                                  'WT_mean': [mean_WT],
                                                  'KO_mean': [mean_KO]}))
    nameStruct_filename = "".join([c for c in nameStruct if c.isalpha() or c.isdigit() or c == ' ']).rstrip()
    mouse_table_nameStruct.to_csv(os.path.join(analysis_path, 'perstructure', nameStruct_filename+'_volumes.csv'))
mouse_table_pername = pd.concat(mouse_table_pername_list, ignore_index=True)
mouse_table_pername = mouse_table_pername.sort_values(by='pVal')
mouse_table_pername.to_csv(os.path.join(analysis_path, 'pername'+'_volumes.csv'))

# Add id_custom column to pVal table
mouse_table_pername = pd.merge(left=mouse_table_pername, right=structure.loc[:, ['name', 'id_custom']],
                               left_on='name', right_on='name')
mouse_table_pername['pVal_inv'] = np.abs(np.log10(mouse_table_pername['pVal']))

# Save separate pval table with only cerebellum
pername_table_cerebellum = pd.merge(left=mouse_table_pername, right=reference_table[['id_custom', 'in_cerebellum']],
         left_on='id_custom', right_on='id_custom')
pername_table_cerebellum = pername_table_cerebellum[pername_table_cerebellum['in_cerebellum']]
pername_table_cerebellum = pername_table_cerebellum.drop('in_cerebellum', 1)
pername_table_cerebellum.to_csv(os.path.join(analysis_path, 'pername_cerebellum_volumes.csv'))



# Create reference images with p-values in the image instead of structure integers
mean_diff = mouse_table_pername['KO_mean'] - mouse_table_pername['WT_mean']
map_from = np.array(mouse_table_pername['id_custom'])
annotation_image = nib.load(annotation_path)
annotation = annotation_image.get_fdata()
for i in range(12):

    annotation_pVal_path = os.path.join(analysis_path, annotation_path.split(os.sep)[-1].split('.')[0])
    if i == 0:
        map_to = np.array(mouse_table_pername['pVal'])
        annotation_pVal_path = annotation_pVal_path + '_pVal' + '.nii.gz'
    elif i == 1:
        map_to = np.array(mouse_table_pername['pVal_inv'])
        annotation_pVal_path = annotation_pVal_path + '_pVal_inv' + '.nii.gz'
    elif i == 2:
        map_to = np.array(mouse_table_pername['pVal']*(mouse_table_pername['pVal'] < 0.05))
        annotation_pVal_path = annotation_pVal_path + '_pVal_sig' + '.nii.gz'
    elif i == 3:
        map_to = np.array(mouse_table_pername['pVal_inv']*(mouse_table_pername['pVal'] < 0.05))
        annotation_pVal_path = annotation_pVal_path + '_pVal_inv_sig' + '.nii.gz'
    elif i == 4:
        map_to = np.array(mouse_table_pername['pVal']*(mean_diff > 0))
        annotation_pVal_path = annotation_pVal_path + '_pVal_volIncrease' + '.nii.gz'
    elif i == 5:
        map_to = np.array(mouse_table_pername['pVal']*(mean_diff < 0))
        annotation_pVal_path = annotation_pVal_path + '_pVal_volDecrease' + '.nii.gz'
    elif i == 6:
        map_to = np.array(mouse_table_pername['pVal_inv']*(mean_diff > 0))
        annotation_pVal_path = annotation_pVal_path + '_pVal_inv_volIncrease' + '.nii.gz'
    elif i == 7:
        map_to = np.array(mouse_table_pername['pVal_inv']*(mean_diff < 0))
        annotation_pVal_path = annotation_pVal_path + '_pVal_inv_volDecrease' + '.nii.gz'
    elif i == 8:
        map_to = np.array(mouse_table_pername['pVal']*(mouse_table_pername['pVal'] < 0.05)*(mean_diff > 0))
        annotation_pVal_path = annotation_pVal_path + '_pVal_sig_volIncrease' + '.nii.gz'
    elif i == 9:
        map_to = np.array(mouse_table_pername['pVal']*(mouse_table_pername['pVal'] < 0.05)*(mean_diff < 0))
        annotation_pVal_path = annotation_pVal_path + '_pVal_sig_volDecrease' + '.nii.gz'
    elif i == 10:
        map_to = np.array(mouse_table_pername['pVal_inv']*(mouse_table_pername['pVal'] < 0.05)*(mean_diff > 0))
        annotation_pVal_path = annotation_pVal_path + '_pVal_inv_sig_volIncrease' + '.nii.gz'
    elif i == 11:
        map_to = np.array(mouse_table_pername['pVal_inv']*(mouse_table_pername['pVal'] < 0.05)*(mean_diff < 0))
        annotation_pVal_path = annotation_pVal_path + '_pVal_inv_sig_volDecrease' + '.nii.gz'

    map_to_filt = np.logical_not(np.isnan(map_to))
    map_to_filtered = map_to[map_to_filt]
    map_from_filtered = map_from[map_to_filt]

    annotation_remapped = np.round(annotation) # always annotation so should never be non-integer
    # input = input.astype(int) # always annotation so should never be non-integer
    annotation_remapped_shape = annotation_remapped.shape
    annotation_remapped = annotation_remapped.reshape(-1)
    annotation_remapped = npi.remap(annotation_remapped, map_from_filtered, map_to_filtered)
    annotation_remapped = annotation_remapped.reshape(annotation_remapped_shape)
    # annotation_remapped = remap_3D(annotation, map_from_filtered.astype(int), map_to_filtered)

    output_image = nib.Nifti1Image(annotation_remapped,
                                   annotation_image.affine)
    nib.save(output_image, annotation_pVal_path)



## Plotting
# cmap = matplotlib.cm.get_cmap('lines')
VOIs = ['Substantia nigra, compact part',
        'Substantia nigra, reticular part',
        'Lobule II',
        'cerebal peduncle',
        'root']
for iVOI in range(len(VOIs)):
    ax = output_table_all[output_table_all['name'] == VOIs[iVOI]][['VolumeNormalized', 'Genotype']].boxplot(by=['Genotype'])

    plt.ylabel('Volume Normalized')
    plt.xlabel('Genotype')
    plt.title(VOIs[iVOI] + ' volumes')
    plt.suptitle('') # that's what you're after
    # ax.set_xticklabels(['WT', 'KO'])
    # plt.show()
    plt.savefig(os.path.join(analysis_path, 'Boxplot_'+VOIs[iVOI]+'_ByGenotype'))




#
#
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




# Create volume tables per structure



# # Load files
# low_detail_cerebellum_image = nib.load(low_detail_cerebellum_path)
# low_detail_cerebellum = low_detail_cerebellum_image.get_fdata()
# high_detail_image = nib.load(high_detail_path)
# high_detail = high_detail_image.get_fdata()
# structure_graph = pd.read_csv(allen_structure_table_path_lowdetail)



# Load remapped structuregraph

# Loop through annotated images

    # Load image

    # Calculate highdetail volumes

    # Infer lowdetail volumes

    # Creat and save table for each mouse

# Combine tables into one great table for all mice in analysis folder

# Plot some boxplot plots






#
# ## Mask volume computation, separate from Allen and invwarped computations
# mouse_voxel_number = mouse_list.copy()
# mouse_mask_volume = mouse_list.copy()
# for iMouse, Mouse in enumerate(mouse_list):
#     mouse_mask_image = nib.load(os.path.join(data_path, Mouse, (Mouse+'_mask_t=500_v=380_k=6.mask.nii.gz')))
#     mouse_mask_image_array = mouse_mask_image.get_fdata()
#
#     mouse_voxel_number[iMouse] = np.sum(mouse_mask_image_array>0)
#     mouse_mask_volume[iMouse] = mouse_voxel_number[iMouse]*voxel_volume
#
# mouse_mask_table = pd.DataFrame({'Mouse': mouse_list, 'MaskVoxelNumber': mouse_voxel_number, 'MaskVolume': mouse_mask_volume})
# mouse_mask_table['Genotype'] = mouse_mask_table['Mouse'].str.split("_", n = 1, expand = True).iloc[:, 0]
# mouse_mask_table['Sex'] = mouse_mask_table['Mouse'].str.split("_", n = 3, expand = True).iloc[:, 2]
# mouse_mask_table.loc[mouse_mask_table['Sex'] != 'female', 'Sex'] = 'male'
#
# # Plotting by genotype
# gen_fig = plt.figure(1)
# annotation_bin = nib.load(allen_image)
# annotation_bin_array = annotation_bin.get_fdata()
# allen_voxelnumber = np.sum(annotation_bin_array>.1)
# allen_volume = allen_voxelnumber * voxel_reference_volume
# mouse_mask_table_plot = pd.concat([mouse_mask_table, pd.DataFrame({'Mouse': ['Allen'], 'MaskVoxelNumber': [allen_voxelnumber], 'MaskVolume': [allen_volume], 'Genotype': ['Allen']})], sort=True)
# # plt.suptitle('figure title')
# # mouse_mask_table_plot=mouse_mask_table
# ax = mouse_mask_table_plot[['MaskVolume', 'Genotype']].boxplot(by=['Genotype'])
# # mouse_mask_table_plot.loc[mouse_mask_table_plot['Genotype'] == 'WT', 'Genotype'] = 1
# # mouse_mask_table_plot.loc[mouse_mask_table_plot['Genotype'] == 'KO', 'Genotype'] = 2
# plt.ylabel('$mm^3$')
# plt.xlabel('Genotype')
# plt.title('Brain mask volumes')
# plt.suptitle('') # that's what you're after
# # ax.set_xticklabels(['WT', 'KO'])
# plt.show()
# plt.savefig(os.path.join(analysis_path, 'Boxplot_MaskVolumes_ByGenotype_Allen'))
#
# # Plotting by genotype and sex
# gensex_fig = plt.figure(2)
# mouse_mask_table[['MaskVolume', 'Genotype', 'Sex']].boxplot(by=['Genotype', 'Sex'])
# plt.ylabel('$mm^3$')
# plt.ylabel('')
# plt.savefig(os.path.join(analysis_path, 'Boxplot_MaskVolumes_ByGenotypeSex'))
# # plt.show()
#
# # pval calculation, equal_var for now ;)
# mouse_mask_table[mouse_mask_table['Genotype'] == 'WT']['MaskVolume']
# cat2 = mouse_mask_table[mouse_mask_table['Genotype'] == 'KO']
#
# print(ttest_ind(mouse_mask_table[mouse_mask_table['Genotype'] == 'WT']['MaskVolume'],
#                 mouse_mask_table[mouse_mask_table['Genotype'] == 'KO']['MaskVolume'],
#                 equal_var=True))
#
# mouse_mask_table.to_csv(os.path.join(analysis_path, 'Mouse_maskvolume_table.csv'))
#
#
#
#
#
#
#
# # Define structure graph
# structure_graph = pd.read_csv(os.path.join(allen_path, 'structure_graph_remapped.csv'))
#
#
#
# ## Get list of invwarped annotation files and compute volumes for them, compiling eventually into single table
# mouse_invwarped_list = glob.glob(data_path+'/*/*/*invwarped*')
# mouse_table_invwarped_list = mouse_invwarped_list.copy()
# for iMouse, Mouse in enumerate(mouse_invwarped_list):
#     Mouse_name = Mouse.split('/')[7]
#     FNIRT_run = Mouse.split('/')[-2]
#
#     mouse_table_invwarped_list[iMouse] = image2volumetable(Mouse, voxel_volume)
#     mouse_table_invwarped_list[iMouse]['Mouse'] = Mouse_name
#     mouse_table_invwarped_list[iMouse]['FNIRT_run'] = FNIRT_run
#
# mouse_table_invwarped = pd.concat(mouse_table_invwarped_list)
#
# mouse_mask_table.to_csv(os.path.join(analysis_path, 'Mouse_invwarped_table.csv'))
#
#
#
# ## Compute volumes for allen reference
# allen_table = image2volumetable(allen_image, voxel_reference_volume)
# allen_table['Mouse'] = 'allen'
# allen_table.to_csv(os.path.join(analysis_path, 'allen_table.csv'))



# ## Compute volumes for flirted allen reference
# allen_flirted_table = image2volumetable(allen_image_flirted, voxel_AMBMC_volume)
# allen_flirted_table['Mouse'] = 'allen_flirted'
# allen_flirted_table.to_csv(os.path.join(analysis_path, 'allen_flirted_table.csv'))
#
#
#
#
#
#
# # for i in range(100):
# #     print(i)
# #     print(np.any(pd.DataFrame(oapi.get_structures(i+1))['id']==182305696))

