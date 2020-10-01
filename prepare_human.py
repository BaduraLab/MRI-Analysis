# make sure qform is the same as sform
# then use fslreorient2std and save archive to Processed

# note that the manual cerebellar adjustment annotation as well as the orsuit annotation need to be downloaded/copied
# manually to each subject folder, automatizing this is probably not worth it.

import nibabel as nib
from glob import glob
import os
from pathlib import Path

# Define paths
raw_path = os.path.join('Data', 'Human', 'Raw')
processed_path = os.path.join('Data', 'Human', 'Processed')
# processed_path = os.path.join('Data', 'Human', 'Processed_EMMetric')
# processed_path = os.path.join('Data', 'Human', 'Processed_SSDMetric')
raw_path_list = glob(os.path.join(raw_path, '*.nii'))

for iRawPath, RawPath in enumerate(raw_path_list):
    raw_name = RawPath.split(os.sep)[-1].split('.')[0]
    raw_foldername = raw_name.split('_')[0]
    Path(os.path.join(processed_path, raw_foldername)).mkdir(parents=True, exist_ok=True)
    ProcessedPath = os.path.join(processed_path, raw_foldername, raw_name+'_reoriented.nii.gz')

    print(RawPath)
    print(ProcessedPath)

    raw_image = nib.load(RawPath)
    processed_image = nib.as_closest_canonical(raw_image)
    processed_image.set_qform(processed_image.affine, code=1)
    processed_image.set_sform(processed_image.affine, code=0)
    nib.save(processed_image, ProcessedPath)
    if 'skull' not in ProcessedPath:
        nib.save(processed_image, ProcessedPath.split('.')[0]+'.nii')
