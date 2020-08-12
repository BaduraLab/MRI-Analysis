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
from scipy import spatial

# Define paths
data_path = os.path.join('Data', 'Human', 'Processed')
reference_path = os.path.join('Data', 'Human', 'Reference')
analysis_path = os.path.join('Data', 'Human', 'Analysis')
input_path_list_list = [glob.glob(os.path.join(data_path, '*', '*annotation_suit_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_subcortical_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_CerebrA_thrarg.nii.gz')),
                        glob.glob(os.path.join(data_path, '*', '*annotation_mask_thrarg.nii.gz'))]
structure_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT.txt'),
                       os.path.join(reference_path, 'subcortical', 'subcortical.csv'),
                       os.path.join(reference_path, 'CerebrA', 'CerebrA.csv'),
                       os.path.join(reference_path, 'CerebrA', 'mask.csv')]
manual_path_list = glob.glob(os.path.join(data_path, '*', '*annotation_suitmask_thrarg_manual.nii.gz'))

# Follows
nAnnotation = len(input_path_list_list)
structure_table_list = [pd.read_csv(structure_path_list[0],
                                    delim_whitespace=True,
                                    names=['VolumeInteger', 'name', 'ID']),
                        pd.read_csv(structure_path_list[1],
                                    names=['VolumeInteger', 'acronym', 'name']),
                        pd.read_csv(structure_path_list[2]),
                        pd.read_csv(structure_path_list[3])]

# Create suit mask for manual adjustment
for iPath in input_path_list_list[0]:
    print(iPath)

    # Get mask output path
    automatic_path = iPath.split('suit')[0] + 'orsuit_thrarg.nii.gz'
    manual_path = iPath.split('suit')[0] + 'suitmask_thrarg_manual.nii.gz'
    adjusted_path = iPath.split('suit')[0] + 'orsuit_thrarg_adjusted.nii.gz'
    print(adjusted_path)

    # Load image
    automatic_image = nib.load(automatic_path)
    automatic = automatic_image.get_fdata()
    manual_image = nib.load(manual_path)
    manual = manual_image.get_fdata()
    adjusted = automatic.copy()

    # Grid for lowdetail and highdetail
    X, Y, Z = np.mgrid[0:automatic_image.shape[0]:1, 0:automatic_image.shape[1]:1, 0:automatic_image.shape[2]:1]

    # Logicals
    orsuit_automatic_logical = automatic > 0
    suit_add_logical = np.logical_and(np.logical_not(orsuit_automatic_logical),
                                      manual > 0)  # no automatic annotation, but there is manual annotation - add
    suit_remove_logical = np.logical_and(np.logical_not(manual),
                                         orsuit_automatic_logical)  # no manual annotation, but there is automatic annotation - remove
    orsuit_automatic_specific_logical = np.isin(automatic, [29, 30, 31, 32, 33, 34])  # specific correct orsuit logical
    suit_remove_logical = np.logical_and(suit_remove_logical, np.logical_not(
        orsuit_automatic_specific_logical))  # do not remove correct orsuit annotation

    # Points
    orsuit_automatic_points = np.vstack((X[orsuit_automatic_logical],
                                         Y[orsuit_automatic_logical],
                                         Z[orsuit_automatic_logical])).transpose()  # Get old automatic points
    suit_add_points = np.vstack((X[suit_add_logical],
                                 Y[suit_add_logical],
                                 Z[suit_add_logical])).transpose()  # Get new manual points

    # Tree
    orsuit_automatic_tree = spatial.KDTree(orsuit_automatic_points) # Get old automatic tree

    # Go through new points
    for iAdd in range(suit_add_points.shape[0]):
        add_index = tuple(suit_add_points[iAdd, :])
        closest_annotated_index = tuple(orsuit_automatic_points[orsuit_automatic_tree.query(add_index)[1]])
        adjusted[add_index] = automatic[closest_annotated_index]

    # Remove points which are not specific orsuit or manual
    adjusted[suit_remove_logical] = 0

    # Save adjusted image
    adjusted_image = nib.Nifti1Image(adjusted,
                                     automatic_image.affine,
                                     automatic_image.header)
    nib.save(adjusted_image, adjusted_path)
