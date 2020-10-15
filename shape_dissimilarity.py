from functions import imageFLIRT2defField
import nibabel as nib
from compress_pickle import dump, load
import numpy as np
import os
import glob

# Define paths
data_path = os.path.join('Data', 'Mouse', 'Processed')
reference_path = os.path.join('Data', 'Mouse', 'Reference')
data_path_list = glob.glob(os.path.join(data_path, '*'))
ref_path = os.path.join(reference_path, 'annotation_25_reoriented_flirted.nii.gz')
ref_image = nib.load(ref_path)

for Path in data_path_list:
    subject = Path.split(os.sep)[-1]
    print(subject)

    # Define subject specific paths
    or_path = os.path.join(Path, subject + '.nii.gz')
    flirtRigid_path = os.path.join(Path, subject + '_flirtRigid.mat')
    flirt_path = os.path.join(Path, subject + '_flirt.mat')
    invflirt_path = os.path.join(Path, subject + '_invflirt.mat')
    defField_path = os.path.join(Path, subject + '_flirt_defField.nii.gz')
    defField_magnitude_path = os.path.join(Path, subject + '_flirt_defField_magnitude.nii.gz')

    or_image = nib.load(or_path)

    with open(flirt_path, 'r') as f:
        txt = f.read()
        flirt = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
    print(flirt)
    with open(invflirt_path, 'r') as f:
        txt = f.read()
        invflirt = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
    print(invflirt)

    # Calculate inverse flirt and check whether it matches with read inverse flirt

    #
    defField = imageFLIRT2defField(ref_image, invflirt)
    defField_image = nib.Nifti1Image(defField, ref_image.affine)
    nib.save(defField_image, defField_path)
    defField_magnitude = np.sqrt(np.sum(np.power(defField, 2), axis=3))
    defField_magnitude_image = nib.Nifti1Image(defField_magnitude, ref_image.affine)
    nib.save(defField_magnitude_image, defField_magnitude_path)



    # Load SyN of subject and get inverse field
    defField_path = os.path.join(Path, subject + '_flirted_syn_defField.nii.gz')
    defField_magnitude_path = os.path.join(Path, subject + '_flirted_syn_defField_magnitude.nii.gz')
    syn_path = os.path.join(Path, subject + '_flirted_syn.pickle.gz')
    with open(syn_path, 'rb') as f:
        [mapping, metric, level_iters, sdr] = load(f, compression='gzip')

    defField = mapping.get_backward_field()
    del mapping, metric, level_iters, sdr

    defField_image = nib.Nifti1Image(defField, ref_image.affine)
    nib.save(defField_image, defField_path)
    defField_magnitude = np.sqrt(np.sum(np.power(defField, 2), axis=3))
    defField_magnitude_image = nib.Nifti1Image(defField_magnitude, ref_image.affine)
    nib.save(defField_magnitude_image, defField_magnitude_path)
