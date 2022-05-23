#!/usr/bin/python3
""" Use RATS generated mask to adjust annotation
"""

__author__ = "Enzo Nio"
__version__ = "1.0.0"
__maintainer__ = "Enzo Nio"

# Mask synned annotation by multiplying annotation with mask
# Note that areas in mask but not in synned annotation are not filled in

import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
# import fsl.wrappers
import csv
from scipy.stats import ttest_ind
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.imwarp import DiffeomorphicMap
from dipy.align.metrics import CCMetric
import datetime
# import SimpleITK as sitk
from compress_pickle import dump, load

# Define
data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
mouse_path_list = glob.glob(os.path.join(data_path, '*'))

# Loop through mice
mouse_path_list = [mouse_path_list[-3]]
for iMousePath, MousePath in enumerate(mouse_path_list):
    ## Adjust synned annotation
    # load images
    Mouse = MousePath.split(os.sep)[-1]
    annotation_path = os.path.join(MousePath, 'allen_annotation_invsynned_to_' + Mouse + '.nii.gz')
    annotation_adjusted_path = os.path.join(MousePath, 'allen_annotation_invsynned_to_' + Mouse + '_adjusted.nii.gz')
    mask_path = os.path.join(MousePath, Mouse + '_mask_t=500_v=380_k=6.mask.nii.gz')

    annotation_image = nib.load(annotation_path)
    annotation_image.get_qform()
    annotation = annotation_image.get_fdata()

    mask_image = nib.load(mask_path)
    mask_image.get_qform()
    mask = mask_image.get_fdata()
    mask = mask / np.max(mask)

    # make annotation 0 everywhere outside of mask
    annotation_adjusted = annotation * mask
    annotation_adjusted_image = nib.Nifti1Image(annotation_adjusted, annotation_image.affine, annotation_image.header)
    print(f'Saving {annotation_adjusted_path}')
    nib.save(annotation_adjusted_image, annotation_adjusted_path)
