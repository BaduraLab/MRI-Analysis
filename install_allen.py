#!/usr/bin/python3

import os, itk
import json
import pandas as pd
from allensdk.api.queries.reference_space_api import ReferenceSpaceApi
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree
from allensdk.config.manifest import Manifest

# the annotation download writes a file, so we will need somewhere to put it
allen_dir = '/home/enzo/Desktop/allen'
Manifest.safe_mkdir(allen_dir)

# this is a string which contains the name of the latest ccf version
allen_version = ReferenceSpaceApi.CCF_VERSION_DEFAULT
allen_resolution = 25 # Set resolution in micrometers

# Download data
rsapi = ReferenceSpaceApi()
# allen_template_path = os.path.join(allen_dir, 'template_'+str(allen_resolution)+'.nrrd')
# rsapi.download_template_volume(resolution=allen_resolution,
#                                file_name=allen_template_path)
allen_annotation_path = os.path.join(allen_dir, 'annotation_'+str(allen_resolution)+'.nrrd')
rsapi.download_annotation_volume(ccf_version=allen_version,
                                 resolution=allen_resolution,
                                 file_name=allen_annotation_path)
allen_average_template_path=os.path.join(allen_dir, 'average_template_'+str(allen_resolution)+'.nrrd')
rsapi.download_volumetric_data(data_path='average_template',
                               coordinate_framework='mouse_ccf',
                               voxel_resolution=allen_resolution,
                               file_name='average_template_'+str(allen_resolution)+'.nrrd',
                               save_file_path=allen_average_template_path)

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

# Save structure graph as json
allen_average_template_json_path=os.path.join(allen_dir, 'structure_graph.json')
with open(allen_average_template_json_path, 'w') as output_json:
    json.dump(allen_structure_graph_dict, output_json)
# Also save structure graph as csv
allen_average_template_csv_path=os.path.join(allen_dir, 'structure_graph.csv')
allen_structure_graph_df = pd.DataFrame(allen_structure_graph_dict)
allen_structure_graph_df.to_csv(allen_average_template_csv_path)



# Convert .nrrd files to .nii
# def convert_image(image_input, image_output):
#     image = itk.imread(image_input)
#     itk.imwrite(image, image_output)

# allen_template_path_nii = os.path.splitext(allen_template_path)[0]+'.nii'
# # convert_image(allen_template_path, allen_template_path_nii)
# os.system('/usr/local/slicer/Slicer '
#           '--no-main-window --no-splash --python-script '
#           './mri_convert_slicer.py '+allen_template_path+' '+allen_template_path_nii)

allen_average_template_path_nii = os.path.splitext(allen_average_template_path)[0]+'.nii'
# convert_image(allen_average_template_path, allen_average_template_path_nii)
os.system('/usr/local/slicer/Slicer '
          '--no-main-window --no-splash --python-script '
          './mri_convert_slicer.py '+allen_average_template_path+' '+allen_average_template_path_nii)

allen_annotation_path_nii = os.path.splitext(allen_annotation_path)[0]+'.nii'
# convert_image(allen_annotation_path, allen_annotation_path_nii)
os.system('/usr/local/slicer/Slicer '
          '--no-main-window --no-splash --python-script '
          './mri_convert_slicer.py '+allen_annotation_path+' '+allen_annotation_path_nii)

# Add
import nibabel as nib
img = nib.load(example_filename)


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