import pandas as pd
import os
import nibabel as nib
import nibabel.processing as nib_processing
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.imwarp import DiffeomorphicMap
from dipy.align.metrics import CCMetric
from compress_pickle import dump, load
import numpy as np
from functions import zeroPadImage
from functions import save_image

reference_path = os.path.join('Data', 'Mouse', 'Reference')
allen_template_path = os.path.join(reference_path, 'average_template_25_reoriented.nii.gz')
allen_template_flirtedRigid_path = os.path.join(reference_path, 'average_template_25_reoriented_flirtedRigid.nii.gz')
allen_template_flirted_path = os.path.join(reference_path, 'average_template_25_reoriented_flirted.nii.gz')
allen_template_flirtRigid_path = os.path.join(reference_path, 'average_template_25_reoriented_flirtRigid.mat')
allen_template_flirt_path = os.path.join(reference_path, 'average_template_25_reoriented_flirt.mat')
allen_template_flirted_synned_path = os.path.join(reference_path,
                                                  'average_template_25_reoriented_flirted_synned.nii.gz')

allen_annotation_path = os.path.join(reference_path, 'annotation_25_reoriented.nii.gz')
allen_annotation_flirtedRigid_path = os.path.join(reference_path, 'annotation_25_reoriented_flirtedRigid.nii.gz')
allen_annotation_flirted_path = os.path.join(reference_path, 'annotation_25_reoriented_flirted.nii.gz')
allen_annotation_flirted_synned_path = os.path.join(reference_path, 'annotation_25_reoriented_flirted_synned.nii.gz')

allen_structure_path = os.path.join(reference_path, 'structure_graph_mc.csv')
allen_template_flirted_syn_path = os.path.join(reference_path, 'average_template_25_reoriented_flirted_syn.pickle.gz')
forw_field_path = os.path.join(reference_path, 'average_template_25_reoriented_flirted_syn_forw.nii.gz')
back_field_path = os.path.join(reference_path, 'average_template_25_reoriented_flirted_syn_back.nii.gz')
AMBMC_template_path = os.path.join(reference_path, 'AMBMC_model.nii.gz')
AMBMC_annotation_path_list = [os.path.join(reference_path, 'AMBMC', 'AMBMC-c57bl6-basalganglia-labels-15um.nii.gz'),
                              os.path.join(reference_path, 'AMBMC', 'AMBMC-c57bl6-cerebellum-labels-15um.nii.gz'),
                              os.path.join(reference_path, 'AMBMC', 'AMBMC-c57bl6-cortex-labels-15um.nii.gz'),
                              os.path.join(reference_path, 'AMBMC', 'AMBMC-c57bl6-hippocampus-labels-15um.nii.gz')]
AMBMC_structure_path_list = [os.path.join(reference_path, 'AMBMC', 'AMBMC_basalganglia_labels.xml'),
                             os.path.join(reference_path, 'AMBMC', 'AMBMC_cerebellum_labels.xml'),
                             os.path.join(reference_path, 'AMBMC', 'AMBMC_cortex_labels.xml'),
                             os.path.join(reference_path, 'AMBMC', 'AMBMC_hippocampus_labels.xml')]

print('reorientation')
correct_qform = nib.load(AMBMC_template_path).get_qform()
AMBMC_path_list = list()
for Path in [AMBMC_template_path] + AMBMC_annotation_path_list:
    print(Path)

    input_image = nib.load(Path)
    input_image.set_qform(correct_qform, 1)

    print(nib.aff2axcodes(input_image.affine))
    output_image = nib.as_closest_canonical(input_image)
    print(nib.aff2axcodes(output_image.affine))

    output_path = Path.split('.')[0] + '_reoriented.nii.gz'
    nib.save(output_image, Path.split('.')[0] + '_reoriented.nii.gz')

    output_resampled_path = Path.split('.')[0] + '_25_reoriented.nii.gz'
    voxel_size = [0.025, 0.025, 0.025]

    output_resampled_image = nib_processing.resample_to_output(output_image, voxel_size)
    nib.save(output_resampled_image, output_resampled_path)

    AMBMC_path_list.append(output_resampled_path)
AMBMC_template_path = AMBMC_path_list[0]

# Load
AMBMC_template_image = nib.load(AMBMC_template_path)
print(nib.aff2axcodes(AMBMC_template_image.affine))
AMBMC_template = AMBMC_template_image.get_fdata()

# AMBMC zero padding
AMBMC_template_zeropadded = zeroPadImage(AMBMC_template, AMBMC_template, 0.1)
AMBMC_template_zeropadded_image = nib.Nifti1Image(AMBMC_template_zeropadded, AMBMC_template_image.affine)
AMBMC_template_zeropadded_path = AMBMC_template_path.split('.')[0]+'_zeropadded.nii.gz'
print(AMBMC_template_zeropadded_path)
nib.save(AMBMC_template_zeropadded_image, AMBMC_template_zeropadded_path)

# FLIRT subject to reference
print('FLIRT rigid start')
os.system('flirt -in ' + allen_template_path + ' \
                 -ref ' + AMBMC_template_zeropadded_path + ' \
                 -out ' + allen_template_flirtedRigid_path + ' \
                 -omat ' + allen_template_flirtRigid_path + ' \
                 -dof ' + '6' + ' \
                 -verbose 0')  # FLIRT subject to reference

# FLIRT allen to AMBMC
print('FLIRT affine start')
os.system('flirt -in ' + allen_template_flirtedRigid_path + ' \
                 -ref ' + AMBMC_template_zeropadded_path + ' \
                 -out ' + allen_template_flirted_path + ' \
                 -omat ' + allen_template_flirt_path + ' \
                 -verbose 0')
allen_template_flirted_image = nib.load(allen_template_flirted_path)
allen_template_flirted = allen_template_flirted_image.get_fdata()

# allen template flirted cropping
allen_template_flirted_cropped = zeroPadImage(allen_template_flirted, allen_template_flirted, 0.1)
allen_template_flirted_cropped_image = nib.Nifti1Image(allen_template_flirted_cropped, allen_template_flirted_image.affine)
allen_template_flirted_cropped_path = allen_template_flirted_path.split('.')[0]+'_cropped.nii.gz'
print(allen_template_flirted_cropped_path)
nib.save(allen_template_flirted_cropped_image, allen_template_flirted_cropped_path)

# SyN flirted images to AMBMC
print('SyN')
metric = CCMetric(3)
level_iters = [10, 10, 5, 5, 5]
sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)
mapping = sdr.optimize(static=AMBMC_template_zeropadded,
                       moving=allen_template_flirted,
                       static_grid2world=AMBMC_template_zeropadded_image.get_qform(),
                       moving_grid2world=allen_template_flirted_image.get_qform())
with open(allen_template_flirted_syn_path, 'wb') as f:
    dump([mapping, metric, level_iters, sdr], f, protocol=4, compression='gzip')
# with open(allen_template_flirted_syn_path, 'rb') as f:
#     [mapping, metric, level_iters, sdr] = load(f, compression='gzip')

forw_field = mapping.get_forward_field()
back_field = mapping.get_backward_field()
forw_SS = np.sum(np.power(forw_field, 2))
back_SS = np.sum(np.power(back_field, 2))
dif_SSD = np.sum(np.power(forw_field + back_field, 2))
dif_SSD_norm = dif_SSD / ((forw_SS + back_SS) / 2)
print(f'dif_SSD_norm = {dif_SSD_norm}')

# forw_field_image = nib.Nifti1Image(forw_field,
#                                    allen_template_flirted_image.affine)
# nib.save(forw_field_image, forw_field_path)
# back_field_image = nib.Nifti1Image(back_field,
#                                    allen_template_flirted_image.affine)
# nib.save(back_field_image, back_field_path)

allen_template_flirted_synned = mapping.transform(allen_template_flirted)
mouse_masked_flirted_synned_image = nib.Nifti1Image(allen_template_flirted_synned,
                                                    allen_template_flirted_image.affine,
                                                    allen_template_flirted_image.header)
nib.save(mouse_masked_flirted_synned_image, allen_template_flirted_synned_path)

# FLIRT (rigid) allen annotation to AMBMC
os.system('flirt -in ' + allen_annotation_path + ' \
                 -ref ' + AMBMC_template_zeropadded_path + ' \
                 -out ' + allen_annotation_flirtedRigid_path + ' \
                 -init ' + allen_template_flirtRigid_path + ' \
                 -dof ' + '6' + ' \
                 -applyxfm' + ' \
                 -verbose 0')

# FLIRT (affine) allen annotation to AMBMC
os.system('flirt -in ' + allen_annotation_flirtedRigid_path + ' \
                 -ref ' + AMBMC_template_zeropadded_path + ' \
                 -out ' + allen_annotation_flirted_path + ' \
                 -init ' + allen_template_flirt_path + ' \
                 -applyxfm' + ' \
                 -verbose 0')
allen_annotation_flirted_image = nib.load(allen_annotation_flirted_path)
allen_annotation_flirted = allen_annotation_flirted_image.get_fdata()

# allen template flirted cropping
allen_annotation_flirted_cropped = zeroPadImage(allen_annotation_flirted, allen_template_flirted, 0.1)
allen_annotation_flirted_cropped_image = nib.Nifti1Image(allen_annotation_flirted_cropped, allen_annotation_flirted_image.affine)
allen_annotation_flirted_cropped_path = allen_annotation_flirted_path.split('.')[0]+'_cropped.nii.gz'
print(allen_annotation_flirted_cropped_path)
nib.save(allen_annotation_flirted_cropped_image, allen_annotation_flirted_cropped_path)

# SyN flirted allen annotation to AMBMC
allen_annotation_flirted_synned = mapping.transform(allen_annotation_flirted)
allen_annotation_flirted_synned_image = nib.Nifti1Image(allen_annotation_flirted_synned,
                                                        allen_annotation_flirted_image.affine,
                                                        allen_annotation_flirted_image.header)
nib.save(allen_annotation_flirted_synned_image, allen_annotation_flirted_synned_path)
