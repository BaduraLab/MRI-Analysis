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
dist = DistanceMetric.get_metric('euclidean')
from sklearn.manifold import MDS
embedding = MDS(n_components=1, dissimilarity='precomputed')
from pathlib import Path
import numpy_indexed as npi



# Define
data_path = os.path.join('Data', 'Human', 'Processed')
reference_path = os.path.join('Data', 'Human', 'Reference')
analysis_path = os.path.join('Data', 'Human', 'Analysis')
structure_MDS_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT_MDS.csv'),
                       os.path.join(reference_path, 'subcortical', 'subcortical_MDS.csv'),
                       os.path.join(reference_path, 'CerebrA', 'CerebrA_MDS.csv')]
structure_table_list = [pd.read_csv(structure_MDS_path_list[0]),
                        pd.read_csv(structure_MDS_path_list[1]),
                        pd.read_csv(structure_MDS_path_list[2])]
annotation_path_list = [os.path.join(reference_path,
                                     'suit',
                                     'atlasesSUIT',
                                     'Lobules-SUIT.nii'),
                        os.path.join(reference_path,
                                     'subcortical',
                                     'prob_atlas_bilateral_thrarg_0.4.nii.gz'),
                        os.path.join(reference_path,
                                     'CerebrA',
                                     'mni_icbm152_CerebrA_tal_nlin_sym_09c_reoriented.nii.gz')]
annotation_name_list = ['suit',
                        'subcortical',
                        'CerebrA']
input_path_list = glob.glob(os.path.join(data_path, '*', 'MDS', '*orsuit*MDS*.nii.gz'))



for Path in input_path_list:
    output_path = Path.split('.')[0] + '_RFMDS.nii.gz'
    print(output_path)

    input_image = nib.load(Path)
    input = input_image.get_fdata()

    input = np.round(input).astype(int)  # ensure integer annotation input

    input_shape = input.shape
    input = input.reshape(-1)
    input = npi.remap(input, list(structure_table_list[0]['VolumeInteger_MDS']), list(structure_table_list[0]['VolumeInteger']))
    input = input.reshape(input_shape)
    output_image = nib.Nifti1Image(input, input_image.affine, input_image.header)
    nib.save(output_image, output_path)
