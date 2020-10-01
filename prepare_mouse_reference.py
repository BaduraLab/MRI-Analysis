import pandas as pd
import os
import nibabel as nib
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.imwarp import DiffeomorphicMap
from dipy.align.metrics import CCMetric
import pickle



reference_path = os.path.join('Data', 'Mouse', 'Reference')
allen_template_path = os.path.join(reference_path, 'average_template_25_reoriented.nii.gz')
allen_annotation_path = os.path.join(reference_path, 'annotation_25_reoriented.nii.gz')
allen_structure_path = os.path.join(reference_path, 'structure_graph_mc.csv')
allen_template_flirted_path = os.path.join(reference_path, 'average_template_25_reoriented_flirted.nii.gz')
allen_template_flirt_path = os.path.join(reference_path, 'average_template_25_reoriented_flirt.mat')
allen_template_flirted_synned_path = os.path.join(reference_path, 'average_template_25_reoriented_flirted_synned.nii.gz')
allen_template_flirted_syn_path = os.path.join(reference_path, 'average_template_25_reoriented_flirted_syn.p')
allen_annotation_flirted_path = os.path.join(reference_path, 'annotation_25_reoriented_flirted.nii.gz')
allen_annotation_flirted_synned_path = os.path.join(reference_path, 'annotation_25_reoriented_flirted_synned.nii.gz')
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

    output_path = Path.split('.')[0]+'_reoriented.nii.gz'
    nib.save(output_image, Path.split('.')[0]+'_reoriented.nii.gz')
    AMBMC_path_list.append(Path)
AMBMC_template_path = AMBMC_path_list[0]



# Load
AMBMC_template_image = nib.load(AMBMC_template_path)
print(nib.aff2axcodes(AMBMC_template_image.affine))
AMBMC_template = AMBMC_template_image.get_fdata()



# FLIRT allen to AMBMC
print('FLIRT')
os.system('flirt -in ' + allen_template_path + ' \
                 -ref ' + AMBMC_template_path + ' \
                 -out ' + allen_template_flirted_path + ' \
                 -omat ' + allen_template_flirt_path + ' \
                 -verbose 1')
allen_template_flirted_image = nib.load(allen_template_flirted_path)
allen_template_flirted = allen_template_flirted_image.get_fdata()



# SyN flirted images to AMBMC
print('SyN')
metric = CCMetric(3)
level_iters = [10, 10, 5, 5, 5]
sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)
mapping = sdr.optimize(static=AMBMC_template,
                       moving=allen_template_flirted,
                       static_grid2world=AMBMC_template_image.get_qform(),
                       moving_grid2world=allen_template_flirted_image.get_qform())
with open(allen_template_flirted_syn_path, 'wb') as f:
    pickle.dump([mapping, metric, level_iters, sdr], f)

allen_template_flirted_synned = mapping.transform(allen_template_flirted)
mouse_masked_flirted_synned_image = nib.Nifti1Image(allen_template_flirted_synned,
                                                    allen_template_flirted_image.affine,
                                                    allen_template_flirted_image.header)
nib.save(mouse_masked_flirted_synned_image, allen_template_flirted_synned_path)



# FLIRT allen annotation to AMBMC
os.system('flirt -in ' + allen_annotation_path + ' \
                 -ref ' + AMBMC_template_path + ' \
                 -out ' + allen_annotation_flirted_path + ' \
                 -init ' + allen_template_flirt_path + ' \
                 -applyfxm' + ' \
                 -verbose 1')
allen_annotation_flirted_image = nib.load(allen_annotation_flirted_path)
allen_annotation_flirted = allen_annotation_flirted_image.get_fdata()



# SyN flirted allen annotation to AMBMC
allen_annotation_flirted_synned = mapping.transform(allen_annotation_flirted)
allen_annotation_flirted_synned_image = nib.Nifti1Image(allen_annotation_flirted_synned,
                                                        allen_annotation_flirted_image.affine,
                                                        allen_annotation_flirted_image.header)
nib.save(allen_annotation_flirted_synned_image, allen_annotation_flirted_synned_path)
