#!/usr/bin/python3
"""  Use manually adjusted human lobular annotations to adjust cerebellar adjusted annotations using nearest neighbour interpolation
"""

__author__ = "Enzo Nio"
__version__ = "1.0.0"
__maintainer__ = "Enzo Nio"

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

# Use _1, adjust it with _2
# For points which are clear in _2 assign second nearest neighbour annotation to _1

# Define paths
data_path = os.path.join('Data', 'Human', 'Processed')
reference_path = os.path.join('Data', 'Human', 'Reference')
analysis_path = os.path.join('Data', 'Human', 'Analysis')
input_path_list_list = [glob.glob(os.path.join(data_path, '*', '*annotation_suit_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_subcortical_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_CerebrA_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_mask_thrarg.nii.gz'))]
structure_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT.csv'),
                       os.path.join(reference_path, 'subcortical', 'subcortical.csv'),
                       os.path.join(reference_path, 'CerebrA', 'CerebrA.csv'),
                       os.path.join(reference_path, 'CerebrA', 'mask.csv')]
structure_table_list = [pd.read_csv(structure_path_list[0]),
                        pd.read_csv(structure_path_list[1]),
                        pd.read_csv(structure_path_list[2]),
                        pd.read_csv(structure_path_list[3])]

volIntList = list(np.arange(34)+1)
volIntList_WM = [29, 30, 31, 32, 33, 34]
volIntList_toconserve = [9]+volIntList_WM
volIntList_toadjust = list(set(volIntList) - set(volIntList_toconserve))



# Adjust orsuit with manual suit - lobular
for iPath in input_path_list_list[0]:
    print(iPath)

    # Get adjusted output path
    original_path = iPath.split('suit')[0] + 'orsuit_thrarg.nii.gz'
    automatic_path = iPath.split('suit')[0] + 'orsuit_thrarg_adjusted_1.nii.gz'
    manual_path = iPath.split('suit')[0] + 'orsuit_thrarg_adjusted_2.nii.gz'
    adjusted_path = iPath.split('suit')[0] + 'orsuit_thrarg_adjusted_lobular.nii.gz'
    print(adjusted_path)

    # Load image
    original_image = nib.load(original_path)
    original = original_image.get_fdata()
    automatic_image = nib.load(automatic_path)
    automatic = automatic_image.get_fdata()
    manual_image = nib.load(manual_path)
    manual = manual_image.get_fdata()

    # Adjust manual with original for volInt 9
    original_9 = original == 9
    manual[np.logical_and(np.logical_not(original_9), manual == 9)] = 0

    adjusted = manual.copy()

    # Grid for lowdetail and highdetail
    X, Y, Z = np.mgrid[0:automatic_image.shape[0]:1, 0:automatic_image.shape[1]:1, 0:automatic_image.shape[2]:1]

    # Logicals
    automatic_logical = np.isin(automatic, volIntList_toadjust)
    manual_logical = np.isin(manual, volIntList_toadjust)
    tree_logical = np.logical_and(manual_logical,
                                 automatic_logical)  # no manual annotation, but there is automatic annotation - add
    add_logical = np.logical_and(np.logical_not(manual_logical),
                                 automatic_logical)  # no manual annotation, but there is automatic annotation - add

    # Points
    tree_points = np.vstack((X[tree_logical],
                                  Y[tree_logical],
                                  Z[tree_logical])).transpose()  # Get old automatic points
    add_points = np.vstack((X[add_logical],
                            Y[add_logical],
                            Z[add_logical])).transpose()  # Get new manual points

    # Tree
    tree = spatial.KDTree(tree_points)  # Get old automatic tree

    # Go through new points
    for iAdd in range(add_points.shape[0]):
        add_index = tuple(add_points[iAdd, :])
        manual_original = manual[add_index]
        automatic_original = automatic[add_index]

        notfound = True
        count = 0
        while notfound:
            count = count + 1
            closest_annotated_indices = tree_points[tree.query(add_index, k=count*1000)[1]]
            for i in range(closest_annotated_indices.shape[0]):
                automatic_new = automatic[tuple(closest_annotated_indices[i,:])]
                if automatic_new != automatic_original:
                    notfound = False
                    adjusted[add_index] = automatic_new
                    # print([manual_original, automatic_original, automatic_new])
                    break

    # Remove points which are not marked in cerebellar mask (_1)
    remove_logical = np.logical_and(np.logical_not(automatic_logical), manual_logical)
    adjusted[remove_logical] = 0

    # Add exception for volInt = 9
    adjusted_extra9 = np.logical_and(np.logical_not(original_9), adjusted == 9)
    if np.any(adjusted_extra9):
        print('volInt=9 added, more 9 in adjusted than expected')
        print(np.sum(adjusted_extra9))
        adjusted[adjusted_extra9] = 0
    adjusted[original_9] = 9

    # Save adjusted image
    adjusted = np.round(adjusted).astype(int)
    adjusted_image = nib.Nifti1Image(adjusted,
                                     manual_image.affine,
                                     manual_image.header)
    nib.save(adjusted_image, adjusted_path)
