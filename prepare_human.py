# make sure qform is the same as sform
# then use fslreorient2std and save archive to Processed

import nibabel as nib
from glob import glob
import os
from pathlib import Path

# Define paths
raw_path = os.path.join('Data', 'Human', 'Raw')
processed_path = os.path.join('Data', 'Human', 'Processed')
raw_path_list = glob(os.path.join(raw_path, '*.nii'))

for iRawPath, RawPath in enumerate(raw_path_list):
    raw_name = RawPath.split(os.sep)[-1].split('.')[0]
    Path(os.path.join(processed_path, raw_name)).mkdir(parents=True, exist_ok=True)
    ProcessedPath = os.path.join(processed_path, raw_name, raw_name+'_reoriented.nii.gz')

    raw_image = nib.load(RawPath)
    processed_image = nib.as_closest_canonical(raw_image)
    processed_image.set_qform(processed_image.affine, code=1)
    processed_image.set_sform(processed_image.affine, code=0)
    nib.save(processed_image, ProcessedPath)
