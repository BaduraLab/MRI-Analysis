#!/usr/bin/python3

import os, itk
import json
import pandas as pd
from allensdk.api.queries.reference_space_api import ReferenceSpaceApi
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree
from allensdk.config.manifest import Manifest
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import nibabel as nib
import numpy as np
import numpy_indexed as npi
import glob
from scipy import spatial



subject = 'KO_6'

# Define paths
allen_fsl_dir = '/usr/local/fsl/data/standard/allen_new'
data_path = '/home/enzo/Desktop/Data/Mouse/Processed_New'
allen_structure_table_path_lowdetail = os.path.join(allen_fsl_dir, 'structure_graph_remapped_lowdetail.csv')
low_detail_cerebellum_path_list = glob.glob(data_path + '/*/*adjusted_lowdetail_manual.nii.gz')
high_detail_path_list = glob.glob(data_path + '/*/*adjusted.nii.gz')
structure_graph = pd.read_csv(allen_structure_table_path_lowdetail)





## Initial variables
# Get all CEREBELLAR highdetail annotation integers (maybe use allen atlas function instead?)
CEREBELLAR_integers_lowdetail = [765, 827, 1186] # 505 and 855 present in high_detail
CEREBELLAR_integers_highdetail = list()
for iRow in range(len(structure_graph)):
    structure_graph_custompath = list(map(int, structure_graph.loc[iRow, 'structure_id_path_custom'].strip('][').split(', ')))
    structure_graph_custompath_inCEREBELLAR = np.isin(structure_graph_custompath, CEREBELLAR_integers_lowdetail)
    if any(structure_graph_custompath_inCEREBELLAR):
        first_inCEREBELLAR_index = np.where(structure_graph_custompath_inCEREBELLAR)[0][0]
        CEREBELLAR_integers_highdetail = CEREBELLAR_integers_highdetail + structure_graph_custompath[first_inCEREBELLAR_index:]
CEREBELLAR_integers_highdetail = np.unique(CEREBELLAR_integers_highdetail)

# Get all cerebellum highdetail annotation integers (maybe use allen atlas function instead?)
cerebellum_integers_lowdetail = 1186
cerebellum_integers_highdetail = list()
for iRow in range(len(structure_graph)):
    structure_graph_custompath = list(map(int, structure_graph.loc[iRow, 'structure_id_path_custom'].strip('][').split(', ')))
    structure_graph_custompath_incerebellum = np.isin(structure_graph_custompath, cerebellum_integers_lowdetail)
    if any(structure_graph_custompath_incerebellum):
        first_incerebellum_index = np.where(structure_graph_custompath_incerebellum)[0][0]
        cerebellum_integers_highdetail = cerebellum_integers_highdetail + structure_graph_custompath[first_incerebellum_index:]
cerebellum_integers_highdetail = np.unique(cerebellum_integers_highdetail)







for iList in range(len(low_detail_cerebellum_path_list)):


    low_detail_cerebellum_path = low_detail_cerebellum_path_list[iList]
    high_detail_path = high_detail_path_list[iList]
    high_detail_cerebellum_path = high_detail_path_list[iList].split('.')[0]+'_cerebellum.nii.gz'


    # Load files
    low_detail_cerebellum_image = nib.load(low_detail_cerebellum_path)
    low_detail_cerebellum = low_detail_cerebellum_image.get_fdata()
    high_detail_image = nib.load(high_detail_path)
    high_detail = high_detail_image.get_fdata()
    high_detail_cerebellum = high_detail.copy()





    # Grid for lowdetail and highdetail
    X, Y, Z = np.mgrid[0:high_detail_image.shape[0]:1, 0:high_detail_image.shape[1]:1, 0:high_detail_image.shape[2]:1]

    # Logicals
    low_detail_cerebellum_logical = low_detail_cerebellum == 1186
    low_detail_not_cerebellum_logical = np.logical_not(low_detail_cerebellum_logical)
    high_detail_CEREBELLAR_logical = np.isin(high_detail, CEREBELLAR_integers_highdetail) # cerebellar commissure, arbor vitae and cerebellum annotation integers respectively (:= CEREBELLAR)
    high_detail_nonzero_logical = high_detail != 0
    high_detail_not_CEREBELLAR_logical = np.logical_and(np.logical_not(high_detail_CEREBELLAR_logical), high_detail_nonzero_logical)

    # points
    high_detail_CEREBELLAR_points = np.vstack((X[high_detail_CEREBELLAR_logical], Y[high_detail_CEREBELLAR_logical], Z[high_detail_CEREBELLAR_logical])).transpose()
    high_detail_not_CEREBELLAR_points = np.vstack((X[high_detail_not_CEREBELLAR_logical], Y[high_detail_not_CEREBELLAR_logical], Z[high_detail_not_CEREBELLAR_logical])).transpose()

    # Trees
    high_detail_CEREBELLAR_tree = spatial.KDTree(high_detail_CEREBELLAR_points)
    high_detail_not_CEREBELLAR_tree = spatial.KDTree(high_detail_not_CEREBELLAR_points)




    # point which are marked as cerebellum in lowdetail, but not annotated as CEREBELLAR in highdetail
    add_high_detail_cerebellum_logical = np.logical_and(low_detail_cerebellum_logical, high_detail_not_CEREBELLAR_logical)
    add_high_detail_cerebellum_points = np.vstack((X[add_high_detail_cerebellum_logical], Y[add_high_detail_cerebellum_logical], Z[add_high_detail_cerebellum_logical])).transpose()

    # for each point annotated as cerebellum in lowdetail, but not annotated as CEREBELLAR in highdetail, find nearest point with CEREBELLAR annotation in highdetail and assign it to an adjusted highdetail
    # fills in a highdetail integers for newly annotated lowdetail cerebellum voxels
    for iAdd in range(len(add_high_detail_cerebellum_points)):
        add_index = tuple(add_high_detail_cerebellum_points[iAdd])
        closest_annotated_index = tuple(high_detail_CEREBELLAR_points[high_detail_CEREBELLAR_tree.query(add_index)[1]])
        high_detail_cerebellum[add_index] = high_detail[closest_annotated_index]



    # points which are marked as not cerebellum in lowdetail, but are marked as CEREBELLAR in highdetail
    remove_high_detail_cerebellum_logical = np.logical_and(low_detail_not_cerebellum_logical, high_detail_CEREBELLAR_logical)
    remove_high_detail_cerebellum_points = np.vstack((X[remove_high_detail_cerebellum_logical], Y[remove_high_detail_cerebellum_logical], Z[remove_high_detail_cerebellum_logical])).transpose()

    # for each point annotated as not cerebellum in lowdetail, but annotated as CEREBELLAR in high detail, find nearest point with non-CEREBELLAR annotation in highdetail and assign it to an adjusted highdetail
    # prevents sometimes accidental removal of non-CEREBELLAR annotation while annotating cerebellum in lowdetail
    for iRemove in range(len(remove_high_detail_cerebellum_points)):
        remove_index = tuple(remove_high_detail_cerebellum_points[iRemove])
        closest_annotated_index = tuple(high_detail_not_CEREBELLAR_points[high_detail_not_CEREBELLAR_tree.query(remove_index)[1]])
        high_detail_cerebellum[remove_index] = high_detail[closest_annotated_index]



    # save adjusted high detail CEREBELLAR adjusted
    high_detail_cerebellum_image = nib.Nifti1Image(high_detail_cerebellum, np.eye(4))
    high_detail_cerebellum_image.set_qform(high_detail_image.get_qform(), code=1)
    high_detail_cerebellum_image.set_sform(np.eye(4), code=0)
    nib.save(high_detail_cerebellum_image, high_detail_cerebellum_path)




    ## Alter high detail with low_detail cerebellum
    # If cerebellum in lowdetail, then cerebellum, arbor vitae or cerebellar commissure in high detail, fill in with nearest neigbour CEREBELLAR annotation

    # If not cerebellum in lowdetail, then not cerebellum in highdetail, fill in values with nearest neighbour non-CEREBELLAR annotation