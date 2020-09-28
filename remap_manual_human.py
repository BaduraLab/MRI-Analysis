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
data_path = os.path.join('Data', 'Human', 'Processed')
reference_path = os.path.join('Data', 'Human', 'Reference')
analysis_path = os.path.join('Data', 'Human', 'Analysis')
annotation_path = os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT.nii')
input_path_list = glob.glob(os.path.join(data_path, '*', '*_lobular.nii.gz'))
# mouse_lobular_path_list = glob.glob(os.path.join(data_path, '*', '*invsynned*cerebellum*lobular.nii.gz'))
reference_structure_path = os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT.csv')
reference_structure_mc_path = os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT_mc.csv')



# Load
structure = pd.read_csv(reference_structure_path)

# Convert reference structure table
id_custom_to_id_mc = [25, 1, 23, 12, 19, 13, 17, 2, 11, 22, 10, 20, 5, 9, 3, 30, 31, 7]
id_mc_to_id_custom = [1, 25, 12, 23, 13, 19, 2, 17, 22, 11, 20, 10, 9, 5, 30, 3, 7, 31]
structure['VolumeInteger_mc'] = npi.remap(structure['VolumeInteger'],
                               id_custom_to_id_mc,
                               id_mc_to_id_custom)
structure.to_csv(reference_structure_mc_path)

# Convert reference and data annotation files
for Path in input_path_list + [annotation_path]:
    annotation_image = nib.load(Path)
    annotation = annotation_image.get_fdata()

    annotation = remap_3D(annotation, id_custom_to_id_mc, id_mc_to_id_custom)

    save_image(annotation, annotation_image, Path.split('.')[0]+'_mc.nii.gz')
