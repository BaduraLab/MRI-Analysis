# prepare/download/preprocess human reference files
import os
import glob
import nibabel as nib
import pandas as pd

# Define paths
reference_path = os.path.join('Data', 'Human', 'Reference')

subcortical_path = os.path.join(reference_path, 'subcortical')
subcortical_image_4D_path = os.path.join(subcortical_path, 'CIT168toMNI152_prob_atlas_bilat_1mm.nii.gz')
subcortical_path_list = glob.glob(os.path.join(subcortical_path, '*volume*.nii.gz'))
subcortical_subcortical_path_list_all = glob.glob(os.path.join(subcortical_path, '*.nii.gz'))

CerebrA_path = os.path.join(reference_path, 'CerebrA')
CerebrA_template_path = os.path.join(CerebrA_path, 'mni_icbm152_t1_tal_nlin_sym_09c.nii')
CerebrA_mask_path = os.path.join(CerebrA_path, 'mni_icbm152_t1_tal_nlin_sym_09c_mask.nii')
CerebrA_masked_path = os.path.join(CerebrA_path, 'mni_icbm152_t1_tal_nlin_sym_09c_masked.nii')
CerebrA_structure_path = os.path.join(CerebrA_path, 'CerebrA_LabelDetails.csv')
CerebrA_structure_adjusted_path = os.path.join(CerebrA_path, 'CerebrA.csv')

# ## subcortical
#
# # convert all subcortical reference files from LAS to RAS
# for iPath, Path in enumerate(subcortical_subcortical_path_list_all):
#     image = nib.load(Path)
#     image = nib.as_closest_canonical(image)
#     image.set_qform(image.affine, code=1)
#     image.set_sform(image.affine, code=0)
#     nib.save(image, Path)
#
#     # transform heeeere to mni space
#     # Path
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





## CerebrA
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
CerebrA_structure_adjusted = CerebrA_structure.rename(columns={'Label Name':'name',
                                                               'RH Label':'Label_RH',
                                                               'LH LabelsNotes':'Label_LH'})
CerebrA_structure_adjusted = pd.wide_to_long(CerebrA_structure_adjusted.loc[:, ['name', 'Label_RH', 'Label_LH']],
                                             stubnames='Label',
                                             i='name',
                                             j='Lateralization',
                                             sep='_', suffix='\w+').reset_index()
CerebrA_structure_adjusted['name_Lateralization'] = CerebrA_structure_adjusted['name'] + ' ' +CerebrA_structure_adjusted['Lateralization']
CerebrA_structure_adjusted = CerebrA_structure_adjusted.loc[:, ['name_Lateralization', 'Label']]
CerebrA_structure_adjusted = CerebrA_structure_adjusted.rename(columns={'name_Lateralization':'name',
                                                                        'Label':'VolumeInteger'})
CerebrA_structure_adjusted.to_csv(CerebrA_structure_adjusted_path)
