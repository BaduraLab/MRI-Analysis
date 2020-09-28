import os
import numpy as np
import nibabel as nib
import pandas as pd
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import glob
import csv
from scipy.stats import ttest_ind
from pathlib import Path
import numpy_indexed as npi
from scipy import spatial
from functions import imageAdjustNN

# Use _2, adjust with _1 which has a certain (cerebellum/sn) mask
# Clear values in _2 as opposed to _1 should be filled in in _2 as NN in _1
# From _2 and _1 _lobular is created

# Define paths
data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
reference_path = os.path.join('Data', 'Mouse', 'Reference')
analysis_path = os.path.join('Data', 'Mouse', 'Analysis')

original_path_list = glob.glob(os.path.join(data_path, '*', '*_cerebellum.nii.gz'))
input1_path_list = glob.glob(os.path.join(data_path, '*', '*_cerebellum_1.nii.gz'))
input2_path_list = glob.glob(os.path.join(data_path, '*', '*_cerebellum_2.nii.gz'))
nInput = len(input1_path_list)
if nInput != len(input2_path_list):
    raise ValueError('Manually adjusted files not in order: len(input1_path_list) != len(input2_path_list)')
# input_path_list = glob.glob(os.path.join(data_path, '*'))
# input_path_list = ['Data\\Mouse\\Processed\\WT_50',
#                    'Data\\Mouse\\Processed\\KO_6',
#                    'Data\\Mouse\\Processed\\KO_2b']
structure_path = os.path.join(reference_path, 'structure_graph_mc.csv')
structure_table = pd.read_csv(structure_path)



# Get cerebellum ids
structure_table['in_cerebellum'] = False
for iVolume in range(structure_table.shape[0]):
    if isinstance(structure_table.loc[iVolume, 'structure_id_path'], str):
        structure_table.loc[iVolume, 'in_cerebellum'] = 512 in list(map(int, structure_table.loc[iVolume, 'structure_id_path'].strip('][').split(', ')))
cerebellum_ids = structure_table[structure_table['in_cerebellum']]['id_custom']
cerebellum_ids = np.round(cerebellum_ids[~np.isnan(cerebellum_ids)]).astype(int)
ce_ids = list(cerebellum_ids) + [827] # 827 is arbor vitae which has been changed manually but not included in cerebellum_ids
sn_ids = [54, 268] # compact and reticular respectively
ce_sn_ids = [ce_ids + sn_ids]




# Adjust MDS_1 orsuit with manual suit
for iPath in range(nInput):
    print(iPath)
    input1_path = input1_path_list[iPath]
    input2_path = input2_path_list[iPath]
    print(input1_path)
    print(input2_path)

    # Get adjusted output path
    adjusted_path = input1_path.split('.')[0][:-2:] + '_lobular.nii.gz'
    print(adjusted_path)

    # Load image
    input1_image = nib.load(input1_path)
    input1 = input1_image.get_fdata()
    input1 = np.round(input1).astype(int)
    input2_image = nib.load(input2_path)
    input2 = input2_image.get_fdata()
    input2 = np.round(input2).astype(int)
    adjusted = input2.copy()

    # Grid for lowdetail and highdetail
    X, Y, Z = np.mgrid[0:input1_image.shape[0]:1, 0:input1_image.shape[1]:1, 0:input1_image.shape[2]:1]

    ## Step 1, fill in cleared lobules in _2 with nearest non-cleared lobuled (determined by _1 mask) in _2

    # Logicals
    input1_logical = np.isin(input1, ce_ids)
    input2_logical = np.isin(input2, ce_ids)
    add_logical = np.logical_and(np.logical_not(input2_logical),
                                 input1_logical)  # no manual annotation, but there is automatic annotation - add
    tree_logical = np.logical_and(input2_logical,
                                  input1_logical)  # no manual annotation, but there is automatic annotation - add

    # Points
    add_points = np.vstack((X[add_logical],
                            Y[add_logical],
                            Z[add_logical])).transpose()  # Get new manual points
    tree_points = np.vstack((X[tree_logical],
                                  Y[tree_logical],
                                  Z[tree_logical])).transpose()  # Get old automatic points

    # Tree
    tree = spatial.KDTree(tree_points)  # Get old automatic tree

    # Go through new points
    for iAdd in range(add_points.shape[0]):
        add_index = tuple(add_points[iAdd, :])
        input1_original = input1[add_index]
        input2_original = input2[add_index]

        notfound = True
        count = 0
        while notfound:
            count = count + 1
            closest_annotated_indices = tree_points[tree.query(add_index, k=count*1000)[1]]
            for i in range(closest_annotated_indices.shape[0]):
                newVal = input2[tuple(closest_annotated_indices[i, :])]
                if newVal != input1_original:
                    notfound = False
                    adjusted[add_index] = newVal
                    # print([manual_original, automatic_original, automatic_new])
                    break

    # Step 2, fill cerebellum in with zeros
    adjusted[np.logical_and(np.isin(adjusted, ce_ids), np.logical_not(np.isin(input1, ce_ids)))] = 0

    # Step 3, fill sn with 54
    adjusted[np.logical_and(np.isin(input1, sn_ids), np.logical_not(np.isin(adjusted, sn_ids)))] = 54

    ## Step 4, fill non-sn with non-sn nearest neighbour

    # Logicals
    input1_sn_logical = np.isin(input1, sn_ids)
    input2_sn_logical = np.isin(input2, sn_ids)
    adjusted_sn_logical = np.isin(adjusted, sn_ids)
    input1_sn_not_logical = np.logical_not(input1_sn_logical)
    input2_sn_not_logical = np.logical_not(input2_sn_logical)
    adjusted_sn_not_logical = np.logical_not(np.isin(adjusted, sn_ids))
    remove_logical = np.logical_and(adjusted_sn_logical, input1_sn_not_logical)
    tree_logical = np.logical_and(input2_sn_not_logical, input1_sn_not_logical)

    # Points
    remove_points = np.vstack((X[remove_logical],
                            Y[remove_logical],
                            Z[remove_logical])).transpose()  # Get new manual points
    tree_points = np.vstack((X[tree_logical],
                             Y[tree_logical],
                             Z[tree_logical])).transpose()  # Get old automatic points

    # Tree
    tree = spatial.KDTree(tree_points)  # Get old automatic tree

    # Go through new points
    for iAdd in range(remove_points.shape[0]):
        add_index = tuple(remove_points[iAdd, :])
        closest_annotated_index = tuple(tree_points[tree.query(add_index)[1]])
        adjusted[add_index] = input2[closest_annotated_index]



    # Save adjusted image
    adjusted_image = nib.Nifti1Image(adjusted,
                                     input2_image.affine,
                                     input2_image.header)
    nib.save(adjusted_image, adjusted_path)
