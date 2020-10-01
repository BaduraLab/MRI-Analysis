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
manual_path_list = glob.glob(os.path.join(data_path, '*', '*annotation_suit_maxprob_thresholded_manual.nii.gz'))
CerebrA_path = os.path.join(reference_path, 'CerebrA')
CerebrA_annotation_path = os.path.join(CerebrA_path, 'mni_icbm152_CerebrA_tal_nlin_sym_09c_reoriented.nii.gz')
CerebrA_cerebellum_volumeIntegers = np.array([46, 97, 2, 53, 20, 71, 50, 101])

# Follows
nAnnotation = len(input_path_list_list)
structure_table_list = [pd.read_csv(structure_path_list[0],
                                    delim_whitespace=True,
                                    names=['VolumeInteger', 'name', 'ID']),
                        pd.read_csv(structure_path_list[1],
                                    names=['VolumeInteger', 'acronym', 'name']),
                        pd.read_csv(structure_path_list[2]),
                        pd.read_csv(structure_path_list[3])]

# Adjust automatic orsuit with manual suit
for iPath in input_path_list_list[0]:
    print(iPath)

    # Get adjusted output path
    automatic_path = iPath.split('suit')[0] + 'orsuit_thrarg.nii.gz'
    manual_path = iPath.split('suit')[0] + 'suit_maxprob_thresholded_manual.nii.gz'
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
    manual_logical = manual > 0
    suit_add_logical = np.logical_and(np.logical_not(orsuit_automatic_logical),
                                      manual_logical)  # no automatic annotation, but there is manual annotation - add
    suit_remove_logical = np.logical_and(np.logical_not(manual_logical),
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
    orsuit_automatic_tree = spatial.KDTree(orsuit_automatic_points)  # Get old automatic tree

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

# Adjust automatic CerebrA with manual suit
for iPath in input_path_list_list[2]:
    print(iPath)

    # Get subject name
    subject_name = iPath.split(os.sep)[-2]

    # Get adjusted output path
    automatic_path = iPath.split('CerebrA')[0] + 'CerebrA_thrarg.nii.gz'
    manual_path = iPath.split('CerebrA')[0] + 'suit_maxprob_thresholded_manual.nii.gz'
    adjusted_path = iPath.split('CerebrA')[0] + 'CerebrA_thrarg_adjusted.nii.gz'
    orsuit_path = iPath.split('CerebrA')[0] + 'orsuit_thrarg.nii.gz'
    print(adjusted_path)

    # Load image
    automatic_image = nib.load(automatic_path)
    automatic = automatic_image.get_fdata()
    automatic = np.squeeze(automatic)
    manual_image = nib.load(manual_path)
    manual = manual_image.get_fdata()
    adjusted = automatic.copy()
    orsuit_image = nib.load(orsuit_path)
    orsuit = orsuit_image.get_fdata()

    # Grid for lowdetail and highdetail
    X, Y, Z = np.mgrid[0:automatic_image.shape[0]:1, 0:automatic_image.shape[1]:1, 0:automatic_image.shape[2]:1]

    # Logicals
    automatic_logical = np.isin(automatic, CerebrA_cerebellum_volumeIntegers)
    manual_logical = manual > 0
    add_logical = np.logical_and(np.logical_not(automatic_logical),
                                 manual_logical)  # no automatic annotation, but there is manual annotation - add
    if subject_name != 'control2':
        remove_logical = np.logical_and(np.logical_not(manual_logical),
                                        automatic_logical)  # no manual annotation, but there is automatic annotation - remove
    else:
        remove_logical = np.logical_and(np.logical_not(manual_logical),
                                        np.isin(adjusted, np.concatenate([CerebrA_cerebellum_volumeIntegers,
                                                                          np.array([39,
                                                                                    90])])))  # no manual annotation, but there is automatic annotation - remove
        orsuit_specific_logical = np.isin(orsuit, [29, 30, 31, 32, 33, 34])
        np.sum(automatic - automatic)
        remove_logical = np.logical_and(remove_logical, np.logical_not(
            orsuit_specific_logical))  # do not remove correct orsuit annotation
    # not_remove_logical = np.logical_not(remove_logical)

    # Points
    automatic_points = np.vstack((X[automatic_logical],
                                  Y[automatic_logical],
                                  Z[automatic_logical])).transpose()  # Get old automatic points
    add_points = np.vstack((X[add_logical],
                            Y[add_logical],
                            Z[add_logical])).transpose()  # Get new manual points
    remove_points = np.vstack((X[remove_logical],
                               Y[remove_logical],
                               Z[remove_logical])).transpose()  # Get remove points

    # Tree
    automatic_tree = spatial.KDTree(automatic_points)  # Get old automatic tree

    # Go through new points
    for iAdd in range(add_points.shape[0]):
        add_index = tuple(add_points[iAdd, :])
        closest_annotated_index = tuple(automatic_points[automatic_tree.query(add_index)[1]])
        adjusted[add_index] = automatic[closest_annotated_index]

    # not remove logical, points and tree
    if subject_name != 'control2':
        not_remove_logical = np.logical_not(np.isin(adjusted,
                                                    CerebrA_cerebellum_volumeIntegers))
    else:
        not_remove_logical = np.logical_not(np.isin(adjusted,
                                                    np.concatenate([CerebrA_cerebellum_volumeIntegers,
                                                                    np.array([39, 90])])))
        print('do not include cerebellar white matter in remove annotation interpolation')
    not_remove_points = np.vstack((X[not_remove_logical],
                                   Y[not_remove_logical],
                                   Z[not_remove_logical])).transpose()  # Get not remove points
    not_remove_tree = spatial.KDTree(not_remove_points)  # Get tree to fill in for removed values (should include 0s)
    not_remove_exception_logical = np.logical_and(not_remove_logical, np.logical_not(adjusted == 0))
    not_remove_exception_points = np.vstack((X[not_remove_exception_logical],
                                             Y[not_remove_exception_logical],
                                            Z[not_remove_exception_logical])).transpose()  # Get not remove expception points
    not_remove_exception_tree = spatial.KDTree(
        not_remove_exception_points)  # Get tree to fill in for removed values (should include 0s)

    # # Remove points
    # adjusted[remove_logical] = 0  # If cerebellum white matter is annotated automatically as cerebellum gray matter
    # # but manually annotated as not cerebellum gray matter, this annotation will be removed and show up as a hole.
    # # You could make an exception for this case in an iRemove loop.

    # Go through remove points
    for iRemove in range(remove_points.shape[0]):
        add_index = tuple(remove_points[iRemove, :])
        closest_annotated_index = tuple(not_remove_points[not_remove_tree.query(add_index)[1]])
        nearest_volumeInteger = adjusted[closest_annotated_index]
        if subject_name != 'control2':
            if nearest_volumeInteger == 0:
                closest_annotated_index_exception = tuple(
                    not_remove_exception_points[not_remove_exception_tree.query(add_index)[1]])
                nearest_volumeInteger_exception = adjusted[closest_annotated_index_exception]
                if nearest_volumeInteger_exception in [39, 90]:
                    adjusted[add_index] = nearest_volumeInteger_exception
                else:
                    adjusted[add_index] = nearest_volumeInteger
            else:
                adjusted[add_index] = nearest_volumeInteger
        else:
            adjusted[add_index] = nearest_volumeInteger
        # you could add an exception for cerebellar white matter where you search a second tree without zeros instead,
        # if this returns cerebellar white matter assign that volumeInteger instead of 0
        # (if zero, search nonzero tree, if cerebellar white matter assign respective integer, if not assign zero)
        # update: exception causes cerebellar white matter annotation to be present clearly far away from cerebellum - turn off exception

    # Make sure that what is marked in orsuit as cerebellar white matter is also marked as cbw in adjusted
    if subject_name == 'control2':
        adjusted = imageAdjustNN(input=adjusted, input_logical=orsuit_specific_logical,
                                 correction=adjusted, correction_logical=np.isin(adjusted, [39, 90]))

    # Save adjusted image
    adjusted_image = nib.Nifti1Image(adjusted,
                                     automatic_image.affine,
                                     automatic_image.header)
    nib.save(adjusted_image, adjusted_path)
