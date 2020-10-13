from functions import imageFLIRT2defField
import nibabel as nib
from compress_pickle import dump, load
import numpy as np

or_path = 'Data/Mouse/Processed/WT_50/WT_50.nii.gz'
defField_path = 'Data/Mouse/Processed/WT_50/WT_50_flirt_defField.nii.gz'
defField_magnitude_path = 'Data/Mouse/Processed/WT_50/WT_50_flirt_defField_magnitude.nii.gz'
or_image = nib.load(or_path)
with open('Data/Mouse/Processed/WT_50/WT_50_flirt.mat', 'r') as f:
    txt = f.read()
    flirt = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
    # for num in l:
    #     if len(num)>5:
    #         print(num)
print(flirt)

defField = imageFLIRT2defField(or_image, flirt)
defField_image = nib.Nifti1Image(defField, or_image.affine)
nib.save(defField_image, defField_path)
defField_magnitude = np.sqrt(np.sum(np.power(defField, 2), axis=3))
defField_magnitude_image = nib.Nifti1Image(defField_magnitude, or_image.affine)
nib.save(defField_magnitude_image, defField_magnitude_path)

# Load SyN of WT_50 and get inverse field
ref_path = 'Data/Mouse/Reference/annotation_25_reoriented_flirted_synned.nii.gz'
defField_path = 'Data/Mouse/Processed/WT_50/WT_50_flirted_syn_defField.nii.gz'
defField_magnitude_path = 'Data/Mouse/Processed/WT_50/WT_50_flirted_syn_defField_magnitude.nii.gz'
ref_image = nib.load(ref_path)
with open('Data/Mouse/Processed/WT_50/WT_50_flirted_syn.pickle.gz', 'rb') as f:
    [mapping, metric, level_iters, sdr] = load(f, compression='gzip')

defField = mapping.get_backward_field()
del mapping, metric, level_iters, sdr

defField_image = nib.Nifti1Image(defField, ref_image.affine)
nib.save(defField_image, defField_path)
defField_magnitude = np.sqrt(np.sum(np.power(defField, 2), axis=3))
defField_magnitude_image = nib.Nifti1Image(defField_magnitude, ref_image.affine)
nib.save(defField_magnitude_image, defField_magnitude_path)
