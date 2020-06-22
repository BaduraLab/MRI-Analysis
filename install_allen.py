#!/usr/bin/python3

import os, itk
import json
import pandas as pd
from allensdk.api.queries.reference_space_api import ReferenceSpaceApi
from allensdk.api.queries.image_download_api import ImageDownloadApi
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree
from allensdk.config.manifest import Manifest
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import nibabel as nib
import numpy as np
import numpy_indexed as npi

# the annotation download writes a file, so we will need somewhere to put it
allen_dir = '/home/enzo/Desktop/allen'
allen_fsl_dir = '/usr/local/fsl/data/standard/allen_new'
Manifest.safe_mkdir(allen_dir)

# this is a string which contains the name of the latest ccf version
allen_version = ReferenceSpaceApi.CCF_VERSION_DEFAULT
allen_resolution = 50 # Set resolution in micrometers



# Download data
mcc = MouseConnectivityCache(resolution=allen_resolution)
annot, annot_header = mcc.get_annotation_volume()
template, template_header = mcc.get_template_volume()

# Define paths
allen_annotation_path = os.path.join(allen_dir, 'annotation_'+str(allen_resolution)+'.nii.gz')
allen_annotation_remapped_path = os.path.join(allen_dir, 'annotation_'+str(allen_resolution)+'_remapped.nii.gz')
allen_average_template_path=os.path.join(allen_dir, 'average_template_'+str(allen_resolution)+'.nii.gz')

# Save nifti
qform = np.array([[0, 0, allen_resolution*pow(10, -3), 0], [-allen_resolution*pow(10, -3), 0, 0, 0], [0, -allen_resolution*pow(10, -3), 0, 0], [0, 0, 0, 1]])
img_annotation = nib.Nifti1Image(annot, np.eye(4))
img_average_template = nib.Nifti1Image(template, np.eye(4))
img_annotation.set_qform(qform, code=1)
img_average_template.set_qform(qform, code=1)
img_annotation.set_sform(np.eye(4), code=0)
img_average_template.set_sform(np.eye(4), code=0)
# img_average_template.set_qform(img_average_template_wrongread.get_qform())
nib.save(img_annotation, allen_annotation_path)
nib.save(img_average_template, allen_average_template_path)



# Get structure graph
oapi = OntologiesApi()
allen_structure_graph_dict = oapi.get_structures([1]) # Get structure graph with structure graph id = 1, which is the Mouse Brain Atlas structure graph

# This removes some unused fields returned by the query
allen_structure_graph_dict = StructureTree.clean_structures(allen_structure_graph_dict)

# Get tree
allen_structure_graph_tree = StructureTree(allen_structure_graph_dict)

# now let's take a look at a structure
allen_structure_graph_tree.get_structures_by_name(['Dorsal auditory area'])

# Look at children or parent of structure, important for later (volume calculations)

# Define path of structure graph table
allen_average_template_csv_path=os.path.join(allen_dir, 'structure_graph.csv')
allen_average_template_csv_remapped_path = os.path.join(allen_dir, 'structure_graph_remapped.csv')

# If structure graph already created, simply load old table
if os.path.exists(allen_average_template_csv_remapped_path):
    allen_structure_graph_df = pd.read_csv(allen_average_template_csv_path)
    allen_structure_graph_remapped_df = pd.read_csv(allen_average_template_csv_remapped_path)

# Else create new structure graph
else:
    # Save structure graph as json
    allen_average_template_json_path=os.path.join(allen_dir, 'structure_graph.json')
    with open(allen_average_template_json_path, 'w') as output_json:
        json.dump(allen_structure_graph_dict, output_json)
    # Also save structure graph as csv
    allen_structure_graph_df = pd.DataFrame(allen_structure_graph_dict)
    allen_structure_graph_df.to_csv(allen_average_template_csv_path)

    allen_structure_graph_remapped_df = allen_structure_graph_df.copy()
    allen_structure_graph_remapped_df['id_custom'] = np.random.permutation(len(allen_structure_graph_df))+1
    structure_id_path_custom = list()
    for iRow in range(len(allen_structure_graph_df)):
        structure_id_path_custom.append(list(npi.remap(allen_structure_graph_remapped_df.loc[iRow, 'structure_id_path'],
                                                  list(allen_structure_graph_remapped_df['id']),
                                                  list(allen_structure_graph_remapped_df['id_custom']))))
    allen_structure_graph_remapped_df['structure_id_path_custom'] = structure_id_path_custom
    allen_structure_graph_remapped_df.to_csv(allen_average_template_csv_remapped_path)



# Alter annotation by remapping integers
annot_shape = annot.shape
annot = annot.reshape(-1)
annot = npi.remap(annot, list(allen_structure_graph_remapped_df['id']), list(allen_structure_graph_remapped_df['id_custom']))
annot = annot.reshape(annot_shape)
img_annotation = nib.Nifti1Image(annot, np.eye(4))
img_annotation.set_qform(qform, code=1)
img_annotation.set_sform(np.eye(4), code=0)
nib.save(img_annotation, allen_annotation_remapped_path)







# Convert .nrrd files to .nii
# def convert_image(image_input, image_output):
#     image = itk.imread(image_input)
#     itk.imwrite(image, image_output)

# allen_template_path_nii = os.path.splitext(allen_template_path)[0]+'.nii'
# # convert_image(allen_template_path, allen_template_path_nii)
# os.system('/usr/local/slicer/Slicer '
#           '--no-main-window --no-splash --python-script '
#           './mri_convert_slicer.py '+allen_template_path+' '+allen_template_path_nii)

# allen_average_template_path_nii = os.path.splitext(allen_average_template_path)[0]+'.nii'
# # convert_image(allen_average_template_path, allen_average_template_path_nii)
# os.system('/usr/local/slicer/Slicer '
#           '--no-main-window --no-splash --python-script '
#           './mri_convert_slicer.py '+allen_average_template_path+' '+allen_average_template_path_nii)
#
# allen_annotation_path_nii = os.path.splitext(allen_annotation_path)[0]+'.nii'
# # convert_image(allen_annotation_path, allen_annotation_path_nii)
# os.system('/usr/local/slicer/Slicer '
#           '--no-main-window --no-splash --python-script '
#           './mri_convert_slicer.py '+allen_annotation_path+' '+allen_annotation_path_nii)
#
# # Add
# import nibabel as nib
# img = nib.load(example_filename)


# # Copy orientation information of other reference to these files
# AMBMC_reference='/usr/local/fsl/data/standard/AMBMC_model_reoriented.nii.gz'
# os.system('fslcpgeom '+AMBMC_reference+' '+allen_template_path_nii)
# os.system('fslcpgeom '+AMBMC_reference+' '+allen_average_template_path_nii)
# os.system('fslcpgeom '+AMBMC_reference+' '+allen_annotation_path_nii)
#
# # Correct orientation of .nii files
# os.system('imagename='+allen_template_path_nii+' ; '
#           'fslorient -deleteorient $imagename ; '
#           'fslorient -setqform 0 0 0.0'+str(allen_resolution)+' 0 -0.0'+str(allen_resolution)+' 0 0 0 0 -0.0'+str(allen_resolution)+' 0 0 0 0 0 1 $imagename ; '
#           'fslorient -setqformcode 1 $imagename ; '
#           'fslreorient2std $imagename "${imagename%.*}_reoriented.nii"')
# os.system('imagename='+allen_average_template_path_nii+' ; '
#           'fslorient -deleteorient $imagename ; '
#           'fslorient -setqform 0 0 0.0'+str(allen_resolution)+' 0 -0.0'+str(allen_resolution)+' 0 0 0 0 -0.0'+str(allen_resolution)+' 0 0 0 0 0 1 $imagename ; '
#           'fslorient -setqformcode 1 $imagename ; '
#           'fslreorient2std $imagename "${imagename%.*}_reoriented.nii"')
# os.system('imagename='+allen_annotation_path_nii+' ; '
#           'fslorient -deleteorient $imagename ; '
#           'fslorient -setqform 0 0 0.0'+str(allen_resolution)+' 0 -0.0'+str(allen_resolution)+' 0 0 0 0 -0.0'+str(allen_resolution)+' 0 0 0 0 0 1 $imagename ; '
#           'fslorient -setqformcode 1 $imagename ; '
#           'fslreorient2std $imagename "${imagename%.*}_reoriented.nii"')

#
# # Download data
# image_api = ImageDownloadApi()
# rsapi = ReferenceSpaceApi()
# # allen_template_path = os.path.join(allen_dir, 'template_'+str(allen_resolution)+'.nrrd')
# # rsapi.download_template_volume(resolution=allen_resolution,
# #                                file_name=allen_template_path)
# allen_annotation_path = os.path.join(allen_dir, 'annotation_'+str(allen_resolution)+'.nrrd')
# rsapi.download_annotation_volume(ccf_version='annotation/ccf_2016',
#                                  resolution=allen_resolution,
#                                  file_name=allen_annotation_path)
# allen_average_template_path=os.path.join(allen_dir, 'average_template_'+str(allen_resolution)+'.nrrd')
# rsapi.download_volumetric_data(data_path='average_template',
#                                coordinate_framework='mouse_ccf',
#                                voxel_resolution=allen_resolution,
#                                file_name='average_template_'+str(allen_resolution)+'.nrrd',
#                                save_file_path=allen_average_template_path)



# # load headers of wrongly read files
# allen_annotation__wrongread_path = os.path.join(allen_dir, 'annotation_'+str(allen_resolution)+'_wrongread.nii')
# allen_average_template_wrongread_path=os.path.join(allen_dir, 'template_'+str(allen_resolution)+'_wrongread.nii')
# img_annotation_wrongread = nib.load(allen_annotation__wrongread_path)
# img_average_template_wrongread = nib.load(allen_average_template_wrongread_path)


# # Reorient images to MNI and save as reoriented images
# fslreorient2std_affine = np.array([[0, 0, 1, 0], [-1, 0, 0, 13.175], [0, -1, 0, 7.975], [0, 0, 0, 1]])
# img_annotation.set_qform(qform*fslreorient2std_affine)
# img_average_template.set_qform(qform*fslreorient2std_affine)
# allen_fsl_annotation_path = os.path.join(allen_fsl_dir, 'annotation_'+str(allen_resolution)+'_reoriented.nii.gz')
# allen_fsl_average_template_path=os.path.join(allen_fsl_dir, 'average_template_'+str(allen_resolution)+'_reoriented.nii.gz')
# nib.save(img_annotation, allen_fsl_annotation_path)
# nib.save(img_average_template, allen_fsl_average_template_path)