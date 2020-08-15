# prepare/download/preprocess human reference files
import os
import glob
import nibabel as nib
import pandas as pd
from pathlib import Path
import numpy as np
import nibabel as nib

# Define paths
reference_path = os.path.join('Data', 'Human', 'Reference')
reference_path_list = list(set(Path(reference_path).rglob('*.nii')) |
                           set(Path(reference_path).rglob('*.nii.gz')))

suit_path = os.path.join(reference_path, 'suit', 'atlasesSUIT')
suit_path_list = glob.glob(os.path.join(suit_path, '*.nii.gz'))+glob.glob(os.path.join(suit_path, '*.nii'))
suit_path_list = list(set(suit_path_list) - set(glob.glob(os.path.join(suit_path, '*reoriented.nii.gz'))))
suit_path_list = suit_path_list + [os.path.join(reference_path, 'standard', 'MNI152_T1_1mm_brain.nii.gz'),
                                   os.path.join(reference_path, 'atlases', 'Cerebellum', 'Cerebellum-MNIfnirt-prob-1mm.nii.gz')]

subcortical_path = os.path.join(reference_path, 'subcortical')
subcortical_path_list = glob.glob(os.path.join(subcortical_path, '*.nii.gz'))
subcortical_path_list = list(set(subcortical_path_list) - set(glob.glob(os.path.join(subcortical_path, '*reoriented.nii.gz'))))
subcortical_annotation_4D_path = os.path.join(subcortical_path, 'prob_atlas_bilateral.nii.gz')
probability_threshold = 0.4
subcortical_annotation_path = os.path.join(subcortical_path,
                                           'prob_atlas_bilateral_thrarg_' + str(probability_threshold) + '.nii.gz')
subcortical_maxprob_path = os.path.join(subcortical_path,
                                        'prob_atlas_bilateral_maxprob.nii.gz')
subcortical_path_list = glob.glob(os.path.join(subcortical_path, '*volume*.nii.gz'))
subcortical_subcortical_path_list_all = glob.glob(os.path.join(subcortical_path, '*.nii.gz'))

CerebrA_path = os.path.join(reference_path, 'CerebrA')
CerebrA_path_list = glob.glob(os.path.join(CerebrA_path, '*.nii'))
CerebrA_template_path = os.path.join(CerebrA_path, 'mni_icbm152_t1_tal_nlin_sym_09c.nii')
CerebrA_mask_path = os.path.join(CerebrA_path, 'mni_icbm152_t1_tal_nlin_sym_09c_mask.nii')
CerebrA_masked_path = os.path.join(CerebrA_path, 'mni_icbm152_t1_tal_nlin_sym_09c_masked.nii')
CerebrA_structure_path = os.path.join(CerebrA_path, 'CerebrA_LabelDetails.csv')
CerebrA_structure_adjusted_path = os.path.join(CerebrA_path, 'CerebrA.csv')
CerebrA_annotation_path = os.path.join(CerebrA_path, 'mni_icbm152_CerebrA_tal_nlin_sym_09c_reoriented.nii.gz')
CerebrA_cerebellum_path = os.path.join(CerebrA_path,
                                       'mni_icbm152_CerebrA_tal_nlin_sym_09c_cerebellum_reoriented.nii.gz')
CerebrA_cerebellum_volumeIntegers = np.array([46, 97, 2, 53, 20, 71, 50, 101])



## CerebrA
# Create masked image and adjust table
# load images and table
CerebrA_template_image = nib.load(CerebrA_template_path)
CerebrA_template = CerebrA_template_image.get_fdata()
CerebrA_mask_image = nib.load(CerebrA_mask_path)
CerebrA_mask = CerebrA_mask_image.get_fdata()

# apply mask and save masked template
CerebrA_masked = CerebrA_template * CerebrA_mask
CerebrA_masked_image = nib.Nifti1Image(CerebrA_masked,
                                       CerebrA_template_image.affine,
                                       CerebrA_template_image.header)
nib.save(CerebrA_masked_image, CerebrA_masked_path)

# structure table from wide to long format and save
CerebrA_structure = pd.read_csv(CerebrA_structure_path, index_col=False)
CerebrA_structure_adjusted = CerebrA_structure.rename(columns={'Label Name': 'name',
                                                               'RH Label': 'Label_RH',
                                                               'LH LabelsNotes': 'Label_LH'})
CerebrA_structure_adjusted = pd.wide_to_long(CerebrA_structure_adjusted.loc[:, ['name', 'Label_RH', 'Label_LH']],
                                             stubnames='Label',
                                             i='name',
                                             j='Lateralization',
                                             sep='_', suffix='\w+').reset_index()
CerebrA_structure_adjusted['name_Lateralization'] = CerebrA_structure_adjusted['name'] + ' ' + \
                                                    CerebrA_structure_adjusted['Lateralization']
CerebrA_structure_adjusted = CerebrA_structure_adjusted.loc[:, ['name_Lateralization', 'Label']]
CerebrA_structure_adjusted = CerebrA_structure_adjusted.rename(columns={'name_Lateralization': 'name',
                                                                        'Label': 'VolumeInteger'})
CerebrA_structure_adjusted.to_csv(CerebrA_structure_adjusted_path)


# Create cerebellum annotation image
CerebrA_annotation_image = nib.load(CerebrA_annotation_path)
CerebrA_annotation = CerebrA_annotation_image.get_fdata()
CerebrA_cerebellum = CerebrA_annotation * np.isin(CerebrA_annotation, CerebrA_cerebellum_volumeIntegers)
CerebrA_cerebellum_image = nib.Nifti1Image(CerebrA_cerebellum,
                                           CerebrA_annotation_image.affine,
                                           CerebrA_annotation_image.header)
nib.save(CerebrA_cerebellum_image, CerebrA_cerebellum_path)

## subcortical

# convert subcortical 4D probability annotation to 3D thrarg/volumeInteger annotation
subcortical_annotation_4D_image = nib.load(subcortical_annotation_4D_path)
subcortical_annotation_4D = subcortical_annotation_4D_image.get_fdata()
subcortical_annotation_4D_maxprob = np.max(subcortical_annotation_4D, axis=3)
subcortical_annotation_4D_maxprob = subcortical_annotation_4D_maxprob / \
                                    np.max(subcortical_annotation_4D_maxprob)
subcortical_annotation_4D_argprob = np.argmax(subcortical_annotation_4D, axis=3) + 1
subcortical_annotation_4D_thrprob = subcortical_annotation_4D_maxprob > probability_threshold
subcortical_annotation = subcortical_annotation_4D_argprob \
                         * subcortical_annotation_4D_thrprob
# thrarg
subcortical_annotation_image = nib.Nifti1Image(subcortical_annotation,
                                               subcortical_annotation_4D_image.affine,
                                               subcortical_annotation_4D_image.header)
nib.save(subcortical_annotation_image, subcortical_annotation_path)
# maxprob
subcortical_maxprob_image = nib.Nifti1Image(subcortical_annotation_4D_maxprob,
                                            subcortical_annotation_4D_image.affine,
                                            subcortical_annotation_4D_image.header)
nib.save(subcortical_maxprob_image, subcortical_maxprob_path)



# # convert all subcortical reference files from LAS to RAS
# for iPath, Path in enumerate(subcortical_subcortical_path_list_all):
#     image = nib.load(Path)
#     image = nib.as_closest_canonical(image)
#     image.set_qform(image.affine, code=1)
#     image.set_sform(image.affine, code=0)
#     nib.save(image, Path)
#
#     # Path##########
#
# # create 4D probability atlas image
# subcortical_image_list = list()
# for iPath in range(len(subcortical_path_list)):
#     Path = subcortical_path_list[0].split(')')[0].split('(')[0] + '(volume ' + str(iPath + 1) + ').nii.gz'
#
#     print(Path)
#
#     subcortical_image = nib.load(Path)
#
#     print(subcortical_image.get_fdata().shape)
#
#     subcortical_image_list.append(subcortical_image)
#
# subcortical_image_4D = nib.concat_images(subcortical_image_list)
# print(subcortical_image_4D.get_fdata().shape)
# nib.save(subcortical_image_4D, subcortical_image_4D_path)








# All  .nii or .nii.gz files in subcortical folder
# set sform to unknown (code=0) and qform to affine
for iPath, Path in enumerate(sum([suit_path_list, subcortical_path_list, CerebrA_path_list], [])):
    print(Path)
    print(iPath)

    image = nib.load(Path)
    print('canonical')
    image = nib.as_closest_canonical(image)
    print('q')
    image.set_qform(image.affine, code=1)
    print('s')
    image.set_sform(image.affine, code=0)
    print('save')
    nib.save(image, Path.split('.')[0] + '_reoriented.nii.gz')




