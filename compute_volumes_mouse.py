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
from statsmodels.stats.multitest import multipletests
from pathlib import Path
import numpy_indexed as npi
from functions import remap_3D
from functions import save_image

from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.imwarp import DiffeomorphicMap
from dipy.align.metrics import CCMetric
from compress_pickle import dump, load
import pathlib



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
    print(np.sum(mouse_table_reference['Volume'][mouse_table_reference['VolumeInteger'] != 0]))

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
data_new_path = os.path.join('Data', 'Mouse', 'Processed')
reference_path = os.path.join('Data', 'Mouse', 'Reference')
analysis_path = os.path.join('Data', 'Mouse', 'Analysis')
mouse_path_list = glob.glob(os.path.join(data_path, '*', '*invsynned*cerebellum_lobular.nii.gz'))
comb_str = '' #########################################
reference_structure_path = os.path.join(reference_path, 'structure_graph_plus' + comb_str + '.csv')
voxel_reference_volume = 0.000125
voxel_volume = 0.000125
annotation_path = os.path.join(reference_path, 'annotation_50_reoriented.nii.gz')
annotation_image = nib.load(annotation_path)
annotation = annotation_image.get_fdata()
nIterBootstrap = 10000
volInt_hemispheres = [1101, 354, 511, 219, 1041]

# Follows
nMouse = len(mouse_path_list)
structure = pd.read_csv(reference_structure_path)
Path(os.path.join(analysis_path, 'perstructure')).mkdir(parents=True, exist_ok=True)



# Reference volumes
reference_table = image2volumetable(annotation_path, voxel_reference_volume)

# Fill in reference additional reference volumes explicitly
reference_table.loc[reference_table['name'] == 'cerebellum lobules I-III', 'Volume'] = \
    float(reference_table.loc[reference_table['name'] == 'Lingula (I)', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Lobule II', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Lobule III', 'Volume'])
reference_table.loc[reference_table['name'] == 'cerebellum lobules VI-VII', 'Volume'] = \
    float(reference_table.loc[reference_table['name'] == 'Declive (VI)', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Simple lobule', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Folium-tuber vermis (VII)', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Paramedian lobule', 'Volume'])
reference_table.loc[reference_table['name'] == 'cerebellum lobules VIII-IX', 'Volume'] = \
    float(reference_table.loc[reference_table['name'] == 'Pyramus (VIII)', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Copula pyramidis', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Uvula (IX)', 'Volume'])
reference_table.loc[reference_table['name'] == 'cerebellum vestibulo', 'Volume'] = \
    float(reference_table.loc[reference_table['name'] == 'Nodulus (X)', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Paraflocculus', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Flocculus', 'Volume'])
reference_table.loc[reference_table['name'] == 'cerebellum hemispheres', 'Volume'] = \
    float(reference_table.loc[reference_table['name'] == 'Simple lobule', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Paramedian lobule', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Crus 1', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Crus 2', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Copula pyramidis', 'Volume']) + \
    float(reference_table.loc[reference_table['name'] == 'Flocculus', 'Volume'])
annotation_hemispheres = np.where(np.isin(annotation, volInt_hemispheres))
annotation_brain = np.where(annotation > 0)
mid_LR_image = annotation.shape[0]/2
mid_LR_hemispheres = np.mean(annotation_hemispheres[0])
mid_LR_brain = np.mean(annotation_brain[0])
# print(f'Mid LR of hemispheres = {mid_LR_hemispheres}, '
#       f'mid LR of image is {mid_LR_image}, '
#       f'mid LR of brain is {mid_LR_brain}')
reference_table.loc[reference_table['name'] == 'cerebellum left hemisphere', 'Volume'] = \
    np.sum(annotation_hemispheres[0] <= mid_LR_brain) * voxel_reference_volume
reference_table.loc[reference_table['name'] == 'cerebellum right hemisphere', 'Volume'] = \
    np.sum(annotation_hemispheres[0] > mid_LR_brain) * voxel_reference_volume



reference_table.to_csv(os.path.join(analysis_path, 'reference_volumes_mouse'))

reference_table['in_cerebellum'] = False
for iVolume in range(reference_table.shape[0]):
    if isinstance(reference_table.loc[iVolume, 'structure_id_path_custom'], str):
        reference_table.loc[iVolume, 'in_cerebellum'] = 1186 in list(map(int, reference_table.loc[iVolume, 'structure_id_path_custom'].strip('][').split(', ')))
reference_cerebellum_table = reference_table[reference_table['in_cerebellum']][['name', 'acronym', 'id_custom', 'structure_id_path_custom', 'VoxelNumber', 'Volume']]
reference_cerebellum_table.to_csv(os.path.join(analysis_path, 'reference_volumes_cerebellum_mouse.csv'))


reference_table['in_sn'] = False
for iVolume in range(reference_table.shape[0]):
    if isinstance(reference_table.loc[iVolume, 'structure_id_path_custom'], str):
        reference_table.loc[iVolume, 'in_sn'] = 2001 in list(map(int, reference_table.loc[iVolume, 'structure_id_path_custom'].strip('][').split(', ')))
reference_sn_table = reference_table[reference_table['in_sn']][['name', 'acronym', 'id_custom', 'structure_id_path_custom', 'VoxelNumber', 'Volume']]
reference_sn_table.to_csv(os.path.join(analysis_path, 'reference_volumes_sn_mouse.csv'))



mouse_table_list = list()
for iMouse, Mouse in enumerate(mouse_path_list):
    subject = Mouse.split(os.path.sep)[-2]
    print(Mouse)
    print(subject)
    mouse_table = image2volumetable(Mouse, voxel_volume)
    mouse_table.to_csv(os.path.join(analysis_path, subject+'_volumes_mouse.csv'))

    mouse_table['Mouse'] = subject
    mouse_table['Genotype'] = subject.split('_')[0]
    mouse_table['Sex'] = subject.split('_')[-1]
    mouse_table.loc[mouse_table['Sex'] != 'female', 'Sex'] = 'male'
    mouse_table = mouse_table[['Mouse', 'Genotype', 'Sex', 'name', 'acronym', 'id_custom', 'structure_id_path_custom',
                               'VoxelNumber', 'Volume']]

    mouse_table_list.append(mouse_table)

mouse_table_all = pd.concat(mouse_table_list, ignore_index=True)

reference_table['rVolume'] = reference_table['Volume']
output_table_all = pd.merge(left=mouse_table_all,
                            right=reference_table.loc[:, ['id_custom', 'name', 'rVolume', 'include_in_test']],
                            left_on=['id_custom', 'name'],
                            right_on=['id_custom', 'name'],
                            how='right')


####################################################### are subjects the same for index and assigned?
# Fill in reference additional reference volumes explicitly
output_table_all.loc[output_table_all['name'] == 'cerebellum lobules I-III', 'Volume'] = \
    np.array(output_table_all.loc[output_table_all['name'] == 'Lingula (I)', 'Volume']) + \
    np.array(output_table_all.loc[output_table_all['name'] == 'Lobule II', 'Volume']) + \
    np.array(output_table_all.loc[output_table_all['name'] == 'Lobule III', 'Volume'])
output_table_all.loc[output_table_all['name'] == 'cerebellum lobules VI-VII', 'Volume'] = \
    np.array(output_table_all.loc[output_table_all['name'] == 'Declive (VI)', 'Volume']) + \
    np.array(output_table_all.loc[output_table_all['name'] == 'Simple lobule', 'Volume']) + \
    np.array(output_table_all.loc[output_table_all['name'] == 'Folium-tuber vermis (VII)', 'Volume']) + \
    np.array(output_table_all.loc[output_table_all['name'] == 'Paramedian lobule', 'Volume'])
output_table_all.loc[output_table_all['name'] == 'cerebellum lobules VIII-IX', 'Volume'] = \
    np.array(output_table_all.loc[output_table_all['name'] == 'Pyramus (VIII)', 'Volume']) + \
    np.array(output_table_all.loc[output_table_all['name'] == 'Copula pyramidis', 'Volume']) + \
    np.array(output_table_all.loc[output_table_all['name'] == 'Uvula (IX)', 'Volume'])
output_table_all.loc[output_table_all['name'] == 'cerebellum vestibulo', 'Volume'] = \
    np.array(output_table_all.loc[output_table_all['name'] == 'Nodulus (X)', 'Volume']) + \
    np.array(output_table_all.loc[output_table_all['name'] == 'Paraflocculus', 'Volume']) + \
    np.array(output_table_all.loc[output_table_all['name'] == 'Flocculus', 'Volume'])
output_table_all.loc[output_table_all['name'] == 'cerebellum hemispheres', 'Volume'] = \
    np.array(output_table_all.loc[output_table_all['name'] == 'Simple lobule', 'Volume']) + \
    np.array(output_table_all.loc[output_table_all['name'] == 'Paramedian lobule', 'Volume']) + \
    np.array(output_table_all.loc[output_table_all['name'] == 'Crus 1', 'Volume']) + \
    np.array(output_table_all.loc[output_table_all['name'] == 'Crus 2', 'Volume']) + \
    np.array(output_table_all.loc[output_table_all['name'] == 'Copula pyramidis', 'Volume']) # [1101, 354, 511, 219, 1041]
for iMouse, Mouse in enumerate(mouse_path_list):
    subject = Mouse.split(os.path.sep)[-2]
    print(f'LR hemisphere calculation, subject = {subject}')

    # load corrected annotation in native space
    mouse_image = nib.load(Mouse)
    mouse_image_array = mouse_image.get_fdata()

    # load flirt
    flirt_path = os.path.join(data_path, subject, 'FLIRT', subject + '_to_allen_model_warpaffine.mat')
    with open(flirt_path, 'r') as f:
        txt = f.read()
        flirt = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
    # print(flirt)

    # calculate voxel center and flirted coordinate center
    mouse_voxel_image_center = np.array(list(np.array(mouse_image.shape) / 2) + [1])
    mouse_voxel_brain_center = np.array(list(np.mean(np.where((mouse_image_array > 0)), axis=1)) + [1])
    annotation_hemispheres = np.where(np.isin(mouse_image_array.astype(int), volInt_hemispheres))
    mouse_voxel_hemisphere_center = np.array(list(np.mean(annotation_hemispheres, axis=1)) + [1])

    mouse_flirted_image_center = flirt.dot(mouse_image.affine.dot(mouse_voxel_image_center))
    mouse_flirted_brain_center = flirt.dot(mouse_image.affine.dot(mouse_voxel_brain_center))
    mouse_flirted_hemisphere_center = flirt.dot(mouse_image.affine.dot(mouse_voxel_hemisphere_center))

    # print(f'center_voxel_image = {mouse_voxel_image_center}, ')
    # print(f'center_voxel_brain = {mouse_voxel_brain_center}, ')
    # print(f'center_voxel_hemispheres = {mouse_voxel_hemisphere_center}')
    #
    # print(f'center_flirted_image = {mouse_flirted_image_center}, ')
    # print(f'center_flirted_brain = {mouse_flirted_brain_center}, ')
    # print(f'center_flirted_hemispheres = {mouse_flirted_hemisphere_center}')



    annotation_hemispheres_L = np.empty(len(annotation_hemispheres[0]))
    annotation_hemispheres_L[:] = False
    for iPoint in range(len(annotation_hemispheres[0])):
        iPoint_voxCoord = np.array([annotation_hemispheres[0][iPoint],
                                  annotation_hemispheres[1][iPoint],
                                  annotation_hemispheres[2][iPoint], 1])
        iPoint_Coord = mouse_image.affine.dot(iPoint_voxCoord)
        # iPoint_Coord[3] = 1
        iPoint_flirted_Coord = flirt.dot(iPoint_Coord)

        # print(f'mouse_flirted_brain_center = {mouse_flirted_brain_center}, ')
        # print(f'iPoint_flirted_Coord = {iPoint_flirted_Coord}, ')
        annotation_hemispheres_L[iPoint] = iPoint_flirted_Coord[0] < mouse_flirted_brain_center[0]
    annotation_hemispheres_R = np.logical_not(annotation_hemispheres_L)

    output_table_all.loc[np.logical_and(output_table_all['name'] == 'cerebellum left hemisphere',
                                        output_table_all['Mouse'] == subject), 'Volume'] = \
        np.sum(annotation_hemispheres_L) * voxel_reference_volume
    output_table_all.loc[np.logical_and(output_table_all['name'] == 'cerebellum right hemisphere',
                                        output_table_all['Mouse'] == subject), 'Volume'] = \
        np.sum(annotation_hemispheres_R) * voxel_reference_volume



    # annotation_hemispheres_flirtedRigid = annotation_hemispheres[0]
    # annotation_flirted_path = os.path.join(data_path, subject, 'FLIRT', subject + '_inmasked_flirted.nii.gz') ######## different reference used! change to reference
    # annotation_flirted_image = nib.load(annotation_flirted_path)
    # annotation_flirtedRigid_path = os.path.join(data_new_path, subject, subject + '_annotation_flirtedRigid.nii.gz') ######## different reference used! change to reference
    # annotation_flirtedRigid_image = nib.load(annotation_flirtedRigid_path)
    # annotation_or_path = os.path.join(data_new_path, subject, subject + '_annotation.nii.gz') ######## different reference used! change to reference
    # annotation_or_image = nib.load(annotation_or_path)
    # annotation_flirted = annotation_flirted_image.get_fdata()
    # print(annotation_flirted_image.affine)
    # annotation_flirted_hemispheres = np.where(np.isin(annotation_flirted.astype(int), volInt_hemispheres))
    # LR_center_flirted_image = annotation_flirted_image.affine.dot(
    #     np.array(list(np.array(annotation_flirted_image.shape) / 2) + [1]))
    # LR_center_flirted_brain = annotation_flirted_image.affine.dot(
    #     np.array(list(np.mean(np.where(annotation_flirted > 0), axis=1)) + [1]))
    # LR_center_flirted_hemispheres = annotation_flirted_image.affine.dot(
    #     np.array(list(np.mean(annotation_flirted_hemispheres, axis=1)) + [1]))

    # flirtRigid_path = os.path.join(data_new_path, subject, subject + '_flirtRigid.mat')
    # with open(flirtRigid_path, 'r') as f:
    #     txt = f.read()
    #     flirtRigid = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
    # print(flirtRigid)

    # mouse_masked_flirted_syn_path = os.path.join(data_new_path, subject, subject + '_flirted_syn.pickle.gz')
    # with open(mouse_masked_flirted_syn_path, 'rb') as f:
    #     [mapping, metric, level_iters, sdr] = load(f, compression='gzip')



output_table_all['VolumeNormalized'] = output_table_all['Volume']/output_table_all['rVolume']

# For each subject add mask volumes for  normalization
output_table_all['rootVolume'] = np.nan
for subject in np.unique(output_table_all['Mouse']):
    output_table_all.loc[output_table_all['Mouse'] == subject, 'rootVolume'] = float(output_table_all.loc[np.logical_and(output_table_all['Mouse'] == subject, output_table_all['name'] == 'root'), 'Volume'])

output_table_all.to_csv(os.path.join(analysis_path, 'all_volumes_mouse.csv'))

output_table_all['VolumeRootNormalized'] = output_table_all['Volume']/output_table_all['rootVolume']

output_table_cerebellum = pd.merge(left=output_table_all, right=reference_table[['id_custom', 'in_cerebellum']],
         left_on='id_custom', right_on='id_custom')
output_table_cerebellum = output_table_cerebellum[output_table_cerebellum['in_cerebellum']]
output_table_cerebellum = output_table_cerebellum[output_table_cerebellum['rVolume'] != 0]
output_table_cerebellum = output_table_cerebellum.drop('in_cerebellum', 1)
output_table_cerebellum.to_csv(os.path.join(analysis_path, 'cerebellum_volumes_mouse.csv'))

output_table_cerebellum = pd.merge(left=output_table_all, right=reference_table[['id_custom', 'in_sn']],
         left_on='id_custom', right_on='id_custom')
output_table_cerebellum = output_table_cerebellum[output_table_cerebellum['in_sn']]
output_table_cerebellum = output_table_cerebellum.drop('in_sn', 1)
output_table_cerebellum.to_csv(os.path.join(analysis_path, 'sn_volumes_mouse.csv'))


# In volume table, go through each structure and determine the p-value between genotypes, create a new p-value table
pername_table_path_list = [os.path.join(analysis_path, 'pername_1_volumes_mouse' + comb_str + '.csv'),
                           os.path.join(analysis_path, 'pername_2_volumes_mouse' + comb_str + '.csv')]
pername_cerebellum_table_path_list = [os.path.join(analysis_path, 'pername_1_cerebellum_volumes_mouse' + comb_str + '.csv'),
                                      os.path.join(analysis_path, 'pername_2_cerebellum_volumes_mouse' + comb_str + '.csv')]
pername_sn_table_path_list = [os.path.join(analysis_path, 'pername_1_sn_volumes_mouse' + comb_str + '.csv'),
                              os.path.join(analysis_path, 'pername_2_sn_volumes_mouse' + comb_str + '.csv')]
for iIncludeInTest in [1, 2]:
    mouse_table_pername_list = list()
    mouse_table_all_nobackground = output_table_all.loc[np.logical_not(pd.isnull(mouse_table_all['name']))]
    mouse_table_all_nobackground = mouse_table_all_nobackground[mouse_table_all_nobackground['include_in_test'] == iIncludeInTest]
    # mouse_table_all.loc[pd.isnull(mouse_table_all['name']), 'name'] = 'background'
    name_uniq = np.unique(np.array(mouse_table_all_nobackground['name'].astype('category')))
    nName = len(name_uniq)
    for nameStruct in name_uniq:
        mouse_table_nameStruct = mouse_table_all_nobackground[mouse_table_all_nobackground['name'] == nameStruct]
        mouse_table_nameStruct_WT = mouse_table_nameStruct.loc[mouse_table_all_nobackground['Genotype'] == 'WT']
        mouse_table_nameStruct_KO = mouse_table_nameStruct.loc[mouse_table_all_nobackground['Genotype'] == 'KO']
        [t_stat, p_val] = ttest_ind(mouse_table_nameStruct_WT['Volume'],
            mouse_table_nameStruct_KO['Volume'],
            equal_var=False)
        [t_stat_RN, p_val_RN] = ttest_ind(mouse_table_nameStruct_WT['VolumeRootNormalized'],
            mouse_table_nameStruct_KO['VolumeRootNormalized'],
            equal_var=False)
        mean_WT = np.mean(mouse_table_nameStruct_WT['Volume'])
        mean_KO = np.mean(mouse_table_nameStruct_KO['Volume'])
        std_WT = np.std(mouse_table_nameStruct_WT['Volume'])
        std_KO = np.std(mouse_table_nameStruct_KO['Volume'])
        mean_WT_AN = np.mean(mouse_table_nameStruct_WT['VolumeNormalized'])
        mean_KO_AN = np.mean(mouse_table_nameStruct_KO['VolumeNormalized'])
        std_WT_AN = np.std(mouse_table_nameStruct_WT['VolumeNormalized'])
        std_KO_AN = np.std(mouse_table_nameStruct_KO['VolumeNormalized'])
        mean_WT_RN = np.mean(mouse_table_nameStruct_WT['VolumeRootNormalized'])
        mean_KO_RN = np.mean(mouse_table_nameStruct_KO['VolumeRootNormalized'])
        std_WT_RN = np.std(mouse_table_nameStruct_WT['VolumeRootNormalized'])
        std_KO_RN = np.std(mouse_table_nameStruct_KO['VolumeRootNormalized'])
        nWT = len(mouse_table_nameStruct_WT['Volume'])
        nKO = len(mouse_table_nameStruct_KO['Volume'])
        N = nWT+nKO
        # print(f'nWT{nWT}, nKO={nKO}, N={N}')
        S_p = np.sqrt(((nWT - 1) * np.power(std_WT, 2) + (nKO - 1) * np.power(std_KO, 2))/(N - 2))
        S_p_AN = np.sqrt(((nWT - 1) * np.power(std_WT_AN, 2) + (nKO - 1) * np.power(std_KO_AN, 2))/(N - 2))
        S_p_RN = np.sqrt(((nWT - 1) * np.power(std_WT_RN, 2) + (nKO - 1) * np.power(std_KO_RN, 2))/(N - 2))
        cohenD = (mean_WT - mean_KO) / S_p
        cohenD_AN = (mean_WT_AN - mean_KO_AN) / S_p_AN
        cohenD_RN = (mean_WT_RN - mean_KO_RN) / S_p_RN

        S_a = np.sqrt((np.power(std_WT, 2) + np.power(std_KO, 2)) / 2)
        S_a_RN = np.sqrt((np.power(std_WT_RN, 2) + np.power(std_KO_RN, 2)) / 2)
        cohenD_RN_a = (mean_WT_RN - mean_KO_RN) / S_a_RN
        hedge_correction = (N - 3) / (N - 2.25)
        cohenD_ac = ((mean_WT - mean_KO) / S_a) * hedge_correction
        cohenD_RN_ac = ((mean_WT_RN - mean_KO_RN) / S_a_RN) * hedge_correction
        # print(f'cohenD_RN = {cohenD_RN}, cohenD_RN_a = {cohenD_RN_a}, '
        #       f'cohenD_RN_ac = {cohenD_RN_ac}, cohenD_ac = {cohenD_ac}')

        cohenD_BS = np.empty(nIterBootstrap)
        cohenD_RN_BS = np.empty(nIterBootstrap)
        for iBS in range(nIterBootstrap):
            WT_BS = np.random.normal(mean_WT, std_WT, nWT)
            KO_BS = np.random.normal(mean_KO, std_KO, nKO)
            mean_WT_BS = np.mean(WT_BS)
            mean_KO_BS = np.mean(KO_BS)
            std_WT_BS = np.std(WT_BS)
            std_KO_BS = np.std(KO_BS)
            S_p = np.sqrt(((nWT - 1) * np.power(std_WT_BS, 2) + (nKO - 1) * np.power(std_KO_BS, 2))/(nKO + nWT - 2))
            S_a = np.sqrt((np.power(std_WT, 2) + np.power(std_KO, 2)) / 2)
            cohenD_BS[iBS] = ((mean_WT_BS - mean_KO_BS)/S_a)*hedge_correction

            WT_RN_BS = np.random.normal(mean_WT_RN, std_WT_RN, nWT)
            KO_RN_BS = np.random.normal(mean_KO_RN, std_KO_RN, nKO)
            mean_WT_RN_BS = np.mean(WT_RN_BS)
            mean_KO_RN_BS = np.mean(KO_RN_BS)
            std_WT_RN_BS = np.std(WT_RN_BS)
            std_KO_RN_BS = np.std(KO_RN_BS)
            S_p_RN = np.sqrt(((nWT - 1) * np.power(std_WT_RN_BS, 2) + (nKO - 1) * np.power(std_KO_RN_BS, 2))/(nKO + nWT - 2))
            S_a_RN = np.sqrt((np.power(std_WT_RN, 2) + np.power(std_KO_RN, 2)) / 2)
            cohenD_RN_BS[iBS] = ((mean_WT_RN_BS - mean_KO_RN_BS)/S_a_RN)*hedge_correction
        cohenD_ac_CI = [np.quantile(cohenD_BS, .025), np.quantile(cohenD_BS, .975)]
        cohenD_RN_ac_CI = [np.quantile(cohenD_RN_BS, .025), np.quantile(cohenD_RN_BS, .975)]
        # cohenD_check = [np.quantile(cohenD_BS, .5), np.median(cohenD_BS)]
        # print(cohenD_check)

        mouse_table_pername_list.append(pd.DataFrame({'name': [nameStruct],
                                                      'cohenD': [cohenD_ac],
                                                      'cohenD_CI': [cohenD_ac_CI],
                                                      'cohenD_BrainNormalized': [cohenD_RN_ac],
                                                      'cohenD_BrainNormalized_CI': [cohenD_RN_ac_CI],
                                                      't_stat': [t_stat],
                                                      'pVal': [p_val],
                                                      'pVal_BrainNormalized': [p_val_RN],
                                                      'WT_mean': [mean_WT],
                                                      'WT_std': [std_WT],
                                                      'KO_mean': [mean_KO],
                                                      'KO_std': [std_KO],
                                                      'WT_mean_AllenNormalized': [mean_WT_AN],
                                                      'WT_std_AllenNormalized': [std_WT_AN],
                                                      'KO_mean_AllenNormalized': [mean_KO_AN],
                                                      'KO_std_AllenNormalized': [std_KO_AN],
                                                      'WT_mean_BrainNormalized': [mean_WT_RN],
                                                      'WT_std_BrainNormalized': [std_WT_RN],
                                                      'KO_mean_BrainNormalized': [mean_KO_RN],
                                                      'KO_std_BrainNormalized': [std_KO_RN]}))

        # nameStruct_filename = "".join([c for c in nameStruct if c.isalpha() or c.isdigit() or c == ' ']).rstrip()
        # mouse_table_nameStruct.to_csv(os.path.join(analysis_path, 'perstructure', nameStruct_filename+'_volumes_mouse.csv'))

    mouse_table_pername = pd.concat(mouse_table_pername_list, ignore_index=True)
    mouse_table_pername['pValBon'] = multipletests(mouse_table_pername['pVal'], method = 'bonferroni')[1]
    mouse_table_pername['pValFDR'] = multipletests(mouse_table_pername['pVal'], method = 'fdr_bh')[1]
    mouse_table_pername.loc[mouse_table_pername['name'] != 'root', 'pValBon_BrainNormalized'] = \
        multipletests(mouse_table_pername.loc[mouse_table_pername['name'] != 'root', 'pVal_BrainNormalized'], method = 'bonferroni')[1]
    mouse_table_pername.loc[mouse_table_pername['name'] != 'root', 'pValFDR_BrainNormalized'] = \
        multipletests(mouse_table_pername.loc[mouse_table_pername['name'] != 'root', 'pVal_BrainNormalized'], method = 'fdr_bh')[1]
    mouse_table_pername = mouse_table_pername.sort_values(by='pVal_BrainNormalized')
    mouse_table_pername = mouse_table_pername.reindex(columns = ['name',
                                                                 'cohenD', 'cohenD_BrainNormalized',
                                                                 'cohenD_CI', 'cohenD_BrainNormalized_CI',
                                                                 't_stat', 'pVal', 'pVal_BrainNormalized',
                                                                 'pValBon', 'pValBon_BrainNormalized',
                                                                 'pValFDR', 'pValFDR_BrainNormalized',
                                                                 'WT_mean', 'WT_std',
                                                                 'KO_mean', 'KO_std',
                                                                 'WT_mean_AllenNormalized', 'WT_std_AllenNormalized',
                                                                 'KO_mean_AllenNormalized', 'KO_std_AllenNormalized',
                                                                 'WT_mean_BrainNormalized', 'WT_std_BrainNormalized',
                                                                 'KO_mean_BrainNormalized', 'KO_std_BrainNormalized'])
    mouse_table_pername.to_csv(pername_table_path_list[iIncludeInTest - 1])

    # Add id_custom column to pVal table
    mouse_table_pername = pd.merge(left=mouse_table_pername, right=structure.loc[:, ['name', 'id_custom']],
                                   left_on='name', right_on='name')
    mouse_table_pername['pVal_inv'] = np.abs(np.log10(mouse_table_pername['pVal']))

    # Save separate pval table with only cerebellum or only sn
    pername_table_cerebellum = pd.merge(left=mouse_table_pername, right=reference_table[['id_custom', 'in_cerebellum', 'in_sn']],
             left_on='id_custom', right_on='id_custom')
    pername_table_cerebellum = pername_table_cerebellum[pername_table_cerebellum['in_cerebellum']]
    pername_table_cerebellum = pername_table_cerebellum.drop('in_cerebellum', 1)
    pername_table_cerebellum.to_csv(pername_cerebellum_table_path_list[iIncludeInTest - 1])

    pername_table_cerebellum = pername_table_cerebellum[pername_table_cerebellum['in_sn']]
    pername_table_cerebellum = pername_table_cerebellum.drop('in_sn', 1)
    pername_table_cerebellum.to_csv(pername_sn_table_path_list[iIncludeInTest - 1])



# Create reference images with p-values in the image instead of structure integers
mean_diff = mouse_table_pername['KO_mean'] - mouse_table_pername['WT_mean']
map_from = np.array(mouse_table_pername['id_custom'])
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
        'root',
        'Folium-tuber vermis (VII)',
        'Lobules IV-V',
        'Pyramus (VIII)',
        'Nucleus accumbens',
        'cerebellum lobules I-III',
        'cerebellum lobules VI-VII',
        'cerebellum lobules VIII-IX',
        'cerebellum vestibulo',
        'cerebellum hemispheres',
        'Vermal regions']
for iVOI in range(len(VOIs)):
    ax = output_table_all[output_table_all['name'] == VOIs[iVOI]][['VolumeNormalized', 'Genotype']].boxplot(by=['Genotype'])

    plt.ylabel('Volume Normalized')
    plt.xlabel('Genotype')
    plt.title(VOIs[iVOI] + ' volumes')
    plt.suptitle('') # that's what you're after
    # ax.set_xticklabels(['WT', 'KO'])
    # plt.show()
    plt.savefig(os.path.join(analysis_path, 'Boxplot_'+VOIs[iVOI]+'_ByGenotype'))

    ax = output_table_all[output_table_all['name'] == VOIs[iVOI]][['VolumeRootNormalized', 'Genotype']].boxplot(by=['Genotype'])

    plt.ylabel('Volume Percentage')
    plt.xlabel('Genotype')
    plt.title(VOIs[iVOI] + ' volumes')
    plt.suptitle('') # that's what you're after
    # ax.set_xticklabels(['WT', 'KO'])
    # plt.show()
    plt.savefig(os.path.join(analysis_path, 'Boxplot_'+VOIs[iVOI]+'_ByGenotype_rootNormalized'))




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

