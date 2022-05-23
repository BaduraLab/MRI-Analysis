#!/usr/bin/python3
"""  Use manually adjusted mouse lobular annotations to adjust cerebellar adjusted annotations using nearest neighbour interpolation
"""

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
# Clear values in _2 as opposed to _1 should be filled in _2 as NN in _1
# From _2 and _1 _lobular is created

# ?? _1 contains added annotation (should overwrite if [sn] or [ce] and [other] NN should be used to fill in _2)
# ?? _2 contains removed annotation (if clear, not [sn] or [ce] and fill in with [other])

# ?? what is the function of having two of these volumes? difference between "fill in clear" and "true clear",
# ?? allows additional adjustment of mask, "true clear" in _1 and "fill in clear" in _2

# Define paths
data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
reference_path = os.path.join('Data', 'Mouse', 'Reference')
analysis_path = os.path.join('Data', 'Mouse', 'Analysis')
mouse_path_list = glob.glob(os.path.join(data_path, '*'))

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
ce_sn_ids = ce_ids + sn_ids
vta_id = [1241]



# Adjust MDS_1 orsuit with manual suit
mouse_path_list = [mouse_path_list[-3]]
for iMousePath, MousePath in enumerate(mouse_path_list):
    # Define paths
    Mouse = MousePath.split(os.sep)[-1]
    original_path = os.path.join(MousePath, 'allen_annotation_invsynned_to_' + Mouse + '_adjusted_cerebellum.nii.gz')
    input1_path = os.path.join(MousePath, 'allen_annotation_invsynned_to_' + Mouse + '_adjusted_cerebellum_1.nii.gz')
    input2_path = os.path.join(MousePath, 'allen_annotation_invsynned_to_' + Mouse + '_adjusted_cerebellum_2.nii.gz')
    print(iMousePath)
    print(MousePath)
    print(original_path)
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

    # Step 2, fill true cleared cerebellum in with zeros
    adjusted[np.logical_and(np.isin(adjusted, ce_ids), np.logical_not(np.isin(input1, ce_ids)))] = 0

    # Step 3, fill sn with 54 -- possible exception for WT_30_female, fill in with nearest neighbour SN
    if Mouse == 'WT_30_female':

         # Logicals
        input1_sn_logical = np.isin(input1, sn_ids)
        input2_sn_logical = np.isin(input2, sn_ids)
        adjusted_sn_logical = np.isin(adjusted, sn_ids)
        input1_sn_not_logical = np.logical_not(input1_sn_logical)
        input2_sn_not_logical = np.logical_not(input2_sn_logical)
        adjusted_sn_not_logical = np.logical_not(np.isin(adjusted, sn_ids))
        remove_logical = np.logical_and(np.logical_and(input2_sn_not_logical, input1_sn_logical),
                                        np.logical_not(np.isin(input2, vta_id)))
        tree_logical = input2_sn_logical

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

    else:
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
    print(f'Saving {adjusted_path}')
    nib.save(adjusted_image, adjusted_path)



# # after that SN needs to be corrected properly for WT_30, with good SNc
# #
# # do this by using SN mask which does not change VTA and NN fills SNr and SNc
# #
# # change old2 to have substantia annotation and be clear in rest
# #
# # cerebellum_2 should contain only the newest annotation of SN and VTA - old2 is template,
# # cerebellum_1 should contain the old manual adjustment -- check
# #
# # this might work for correct SN annotation, although SNc might end up bigger which you could have to correct again
#
# old1_path = mouse_path_list[-3]+'/allen_annotation_invsynned_to_WT_30_female_adjusted_cerebellum_1_old.nii.gz'
# old1_image = nib.load(old1_path)
# nib.aff2axcodes(old1_image.affine)
# old1 = old1_image.get_fdata()
# old1 = old1.astype(int)
#
# old2_path = mouse_path_list[-3]+'/allen_annotation_invsynned_to_WT_30_female_adjusted_cerebellum_2_old.nii.gz'
# old2_image = nib.load(old2_path)
# nib.aff2axcodes(old2_image.affine)
# old2 = old2_image.get_fdata()
# old2 = old2.astype(int)
#
# new1_path = mouse_path_list[-3]+'/allen_annotation_invsynned_to_WT_30_female_adjusted_cerebellum_1_new.nii.gz'
# new1_image = nib.load(new1_path)
# nib.aff2axcodes(new1_image.affine)
# new1 = new1_image.get_fdata()
# new1 = new1.astype(int)
#
# new2_path = mouse_path_list[-3]+'/allen_annotation_invsynned_to_WT_30_female_adjusted_cerebellum_2_new.nii.gz'
# new2_image = nib.load(new2_path)
# nib.aff2axcodes(new2_image.affine)
# new2 = new2_image.get_fdata()
# new2 = new2.astype(int)
#
# # # what is not ce in old2 and is ce in old1 should be clear in new2
# # new2[np.logical_and(np.logical_not(np.isin(old2, ce_ids)), np.isin(old1, ce_ids))] = 0
#
# # replicate old ce in new2 and use old1 as current correction
# new2[np.isin(old2, ce_ids)] = old2[np.isin(old2, ce_ids)]
# new2[np.logical_and(np.logical_not(np.isin(old2, ce_ids)), np.isin(new2, ce_ids))] = 0
#
# # # what is ce in old1 should be ce in new1 --
# # new1[np.isin(old1, ce_ids)] = old1[np.isin(old1, ce_ids)]
# # # what is clear in old1 should be clear in new1 -- just extra mask adjustment
# # new1[old1 == 0] = 0
# # # what is sn in old1 should be sn in new1
# # new1[np.isin(old1, sn_ids)] = old1[np.isin(old1, sn_ids)]
# # new1[np.logical_and(np.isin(new1, sn_ids), np.logical_not(old1, sn_ids))] = 0
#
# # # only update old2 to have SN and VTA of new2
# # old2[np.isin(new2, sn_ids)] = new2[np.isin(new2, sn_ids)]
# # old2[np.logical_and(np.logical_not(np.isin(new2, sn_ids)), np.isin(old2, sn_ids))] = 0
# # old2[np.isin(new2, vta_id)] = new2[np.isin(new2, vta_id)]
# # old2[np.logical_and(np.logical_not(np.isin(new2, vta_id)), np.isin(old2, vta_id))] = 0
#
# # what if ce in new2 but not ce in new1?
#
# # new1_new_path = mouse_path_list[-3]+'/allen_annotation_invsynned_to_WT_30_female_adjusted_cerebellum_1.nii.gz'
# # adjusted_image = nib.Nifti1Image(new1,
# #                                  new1_image.affine,
# #                                  new1_image.header)
# # print(f'Saving {new1_new_path}')
# # nib.save(adjusted_image, new1_new_path)
#
# # old2 = old2.astype(int)
# new2_new_path = mouse_path_list[-3]+'/allen_annotation_invsynned_to_WT_30_female_adjusted_cerebellum_2.nii.gz'
# adjusted_image = nib.Nifti1Image(new2,
#                                  new2_image.affine,
#                                  new2_image.header)
# print(f'Saving {new2_new_path}')
# nib.save(adjusted_image, new2_new_path)
