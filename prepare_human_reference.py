# prepare/download/preprocess human reference files
import os
import glob
import nibabel as nib

reference_path = os.path.join('Data', 'Human', 'Reference')
subcortical_path = os.path.join(reference_path,
    'subcortical')
subcortical_image_4D_path = os.path.join(subcortical_path, 'CIT168toMNI152_prob_atlas_bilat_1mm.nii.gz')

subcortical_path_list = glob.glob(os.path.join(subcortical_path,
                       '*volume*.nii.gz'))

subcortical_image_list = list()
for Path in subcortical_path_list:
    subcortical_image = nib.load(Path)
    subcortical_image_list.append(subcortical_image)

subcortical_image_4D = nib.concat_images(subcortical_image_list)
nib.save(subcortical_image_4D, subcortical_image_4D_path)


# transform all existing templates to MNI space?