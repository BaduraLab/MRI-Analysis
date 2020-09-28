import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import glob
import csv
from scipy.stats import ttest_ind
from scipy import ndimage
from sklearn.neighbors import DistanceMetric
from functions import remap_3D, save_image

dist = DistanceMetric.get_metric('euclidean')
from sklearn.manifold import MDS

embedding = MDS(n_components=1, dissimilarity='precomputed')
from pathlib import Path
import numpy_indexed as npi

# Define
data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
reference_path = os.path.join('Data', 'Mouse', 'Reference')
analysis_path = os.path.join('Data', 'Mouse', 'Analysis')
annotation_path = os.path.join(reference_path, 'annotation_50_reoriented.nii.gz')
mouse_path_list = glob.glob(os.path.join(data_path, '*', '*invsynned*cerebellum.nii.gz'))
mouse_lobular_path_list = glob.glob(os.path.join(data_path, '*', '*invsynned*cerebellum*lobular.nii.gz'))
reference_structure_path = os.path.join(reference_path, 'structure_graph_remapped_lowdetail.csv')
reference_structure_mc_path = os.path.join(reference_path, 'structure_graph_mc.csv')
# reference_structure_output_path = os.path.join(reference_path, 'structure_graph_remapped_lowdetail_mc.csv')
analysis_table_path = os.path.join(analysis_path, 'all_volumes.csv')
# VOIs = ['Lobule II', 'Lobules IV-V', 'Substantia nigra, compact part', 'Substantia nigra, reticular part']
# VOIs = ['Substantia nigra, compact part', 'Substantia nigra, reticular part']



# Load
structure = pd.read_csv(reference_structure_path)
volume_table = pd.read_csv(analysis_table_path)

# Convert reference structure table
id_custom_to_id_mc = [762, 500, 821, 300, 977, 350]
id_mc_to_id_custom = [500, 762, 300, 821, 350, 977]
structure['id_mc'] = npi.remap(structure['id_custom'],
                               id_custom_to_id_mc,
                               id_mc_to_id_custom)
structure.to_csv(reference_structure_mc_path)

# Convert reference and data annotation files
for Path in mouse_path_list + [annotation_path] + mouse_lobular_path_list:
    print(Path)

    annotation_image = nib.load(Path)
    annotation = annotation_image.get_fdata()

    annotation = remap_3D(annotation, id_custom_to_id_mc, id_mc_to_id_custom)

    save_image(annotation, annotation_image, Path.split('.')[0]+'_mc.nii.gz')
