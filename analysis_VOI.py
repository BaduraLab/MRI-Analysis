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
from pathlib import Path



# Define
data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
reference_path = os.path.join('Data', 'Mouse', 'Reference')
analysis_path = os.path.join('Data', 'Mouse', 'Analysis')
mouse_path_list = glob.glob(os.path.join(data_path, '*', '*invsynned*cerebellum.nii.gz'))
reference_structure_path = os.path.join(reference_path, 'structure_graph_remapped_lowdetail.csv')
analysis_table_path = os.path.join(analysis_path, 'all_volumes.csv')
VOIs = ['Lobule II', 'Lobules IV-V', 'Substantia nigra, compact part', 'Substantia nigra, reticular part']
VOIs = ['Substantia nigra, compact part', 'Substantia nigra, reticular part']

# Load
structure = pd.read_csv(reference_structure_path)
volume_table = pd.read_csv(analysis_table_path)
volume_table = pd.merge(left=volume_table, right=structure.loc[:, ['name', 'id_custom']],
                         left_on='name', right_on='name')



# Go through mice and compute center of mass of VOIs
for iMousePath, MousePath in enumerate(mouse_path_list):
    mouse_string = MousePath.split(os.sep)[-2]
    volume_table_mouse = volume_table[volume_table['Mouse'] == mouse_string]

    # Load image
    mouse_image = nib.load(MousePath)
    mouse = mouse_image.get_fdata()

    # Loop over VOIs
    for iVOI, VOI in enumerate(VOIs):
        # Get id_custom of VOI
        mouse_VOI_path = os.path.join(analysis_path, mouse_string + '_' + VOI)
        volume_table_mouse_structure = volume_table_mouse[volume_table_mouse['name'] == VOI]
        VOI_id_custom = volume_table_mouse_structure['id_custom'].iloc[0]

        # for each VOI, determine its center of mass in both affine and voxel coordinates
        mouse_VOI_logical = mouse == VOI_id_custom
        VOI_indices = np.array(np.where(mouse_VOI_logical))
        VOI_indices_center = np.mean(VOI_indices, 1)
        VOI_indices_center_affine = mouse_image.affine.dot(np.concatenate((VOI_indices_center, np.array([0]))))

        print(mouse_string)
        print(VOI)
        print(VOI_indices_center)
        print(VOI_indices_center_affine)

        # for each VOI, isolate it and create new nifti files
        mouse_VOI = mouse.copy()
        mouse_VOI[np.logical_not(mouse_VOI_logical)] = 0
        mouse_VOI_image = nib.Nifti1Image(mouse_VOI, mouse_image.affine, mouse_image.header)
        nib.save(mouse_VOI_image, mouse_VOI_path)
        # view in mricrogl outside of python


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