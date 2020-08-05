# Reorient main (already reoriented) mouse data files so that nibabel is able to read them correctly

import nibabel as nib
import glob
import os

# Define
data_path = os.path.join('Data', 'Mouse', 'Processed')
mouse_path_list = glob.glob(os.path.join(data_path, '*'))

# Loop through mice
for iMousePath, MousePath in enumerate(mouse_path_list):

    # Define mouse paths
    mouse_string = MousePath.split(os.sep)[-1]
    mouse_path = os.path.join(MousePath, mouse_string + '.nii.gz')
    mouse_reoriented_path = os.path.join(MousePath, mouse_string + '_reoriented.nii.gz')

    mouse_image = nib.load(mouse_path)
    nib.aff2axcodes(mouse_image.affine)
    mouse_reoriented_image = nib.as_closest_canonical(mouse_image)
    nib.save(mouse_reoriented_image, mouse_reoriented_path)

    # mask_path = os.path.join(MousePath, mouse_string + '_mask_t=500_v=380_k=6.mask.nii.gz')
    # mouse_masked_path = os.path.join(MousePath, mouse_string + '_masked.nii.gz')
    # mouse_masked_translated_path = os.path.join(MousePath, mouse_string + '_translated.nii.gz')
    # mouse_masked_flirted_path = os.path.join(MousePath, mouse_string + '_flirted.nii.gz')
    # mouse_masked_flirt_path = os.path.join(MousePath, mouse_string + '_flirt.mat')
    # mouse_masked_invflirt_path = os.path.join(MousePath, mouse_string + '_invflirt.mat')
    # mouse_masked_flirted_synned_path = os.path.join(MousePath, mouse_string + '_flirted_synned.nii.gz')
    # reference_annotation_invsynned_path = os.path.join(MousePath, mouse_string + '_annotation_flirted.nii.gz')
    # reference_annotation_invsynned_invflirted_path = os.path.join(MousePath, mouse_string + '_annotation.nii.gz')