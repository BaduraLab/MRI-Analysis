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





# Define paths
allen_fsl_dir = '/usr/local/fsl/data/standard/allen_new'
data_path = '/home/enzo/Desktop/Data/Mouse/Processed_New'
invwarped_list = list(set(glob.glob(data_path+'/*/*invsynned*.nii.gz')) - set(glob.glob(data_path+'/*/*invsynned*flirted.nii.gz')) - set(glob.glob(data_path+'/*/*invsynned*lowdetail.nii.gz')))
# allen_annotation_path = os.path.join(allen_fsl_dir, 'annotation_25_to_AMBMC_flirted.nii.gz')
allen_structure_table_path = os.path.join(allen_fsl_dir, 'structure_graph_remapped.csv')
allen_structure_table_path_lowdetail = os.path.join(allen_fsl_dir, 'structure_graph_remapped_lowdetail.csv')



for iInv in range(len(invwarped_list)):

    # Load files
    invwarped_image = nib.load(invwarped_list[iInv])
    invwarped = invwarped_image.get_fdata()
    structure_graph = pd.read_csv(allen_structure_table_path)

    structure_id_low_detail = list()
    for iRow in range(len(structure_graph)):
        structure_graph_path = list(map(int, structure_graph.loc[iRow, 'structure_id_path_custom'].strip('][').split(', ')))
        structure_graph_path_len = len(structure_graph_path)

        # structure is already low detail (level 3 or lower), remap integer to itself (no change in annotation)
        if structure_graph_path_len < 3:
            structure_id_low_detail.append(structure_graph_path[-1])
        # structure has parent with lower detail, remap integer to lower detail annotation volume
        else:
            if structure_graph_path[2]==90:
                structure_id_low_detail.append(structure_graph_path[-1])
            else:
                structure_id_low_detail.append(structure_graph_path[2])

    structure_graph['id_low_detail'] = structure_id_low_detail
    structure_graph.to_csv(allen_structure_table_path_lowdetail)



    # Alter invwarped by remapping integers
    invwarped_shape = invwarped.shape
    invwarped = invwarped.reshape(-1)
    invwarped = npi.remap(invwarped, list(structure_graph['id_custom']), list(structure_graph['id_low_detail']))
    invwarped = invwarped.reshape(invwarped_shape)
    img_invwarped = nib.Nifti1Image(invwarped, np.eye(4))
    img_invwarped.set_qform(invwarped_image.get_qform(), code=1)
    img_invwarped.set_sform(np.eye(4), code=0)
    invwarped_lowdetail_path = invwarped_list[iInv]
    invwarped_lowdetail_path = os.path.splitext(os.path.splitext(invwarped_lowdetail_path)[0])[0]
    invwarped_lowdetail_path = invwarped_lowdetail_path + '_lowdetail.nii.gz'
    nib.save(img_invwarped, invwarped_lowdetail_path)
