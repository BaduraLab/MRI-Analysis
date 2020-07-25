import os
import nibabel as nib
import datetime
import numpy as np
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.imwarp import DiffeomorphicMap
from dipy.align.metrics import CCMetric
import glob



# Define paths
data_path = os.path.join('Data', 'Human', 'Processed')
mouse_path_list = glob.glob(os.path.join(data_path, '*'))
reference_path = os.path.join('Data', 'Human', 'Reference')
# annotation_path = os.path.join('atlases', 'Cerebellum', 'Talairach', 'Talairach-labels-1mm.nii.gz')
annotation_path = os.path.join('atlases', 'Cerebellum', 'Cerebellum-MNIfnirt-prob-1mm.nii.gz')
template_path = os.path.join('standard', 'MNI152_T1_1mm_brain.nii.gz')
input_path_list = glob.glob(os.path.join(data_path, '*', '*_reoriented.nii.gz'))



# Loop through inputs
for iInputPath, InputPath in enumerate(input_path_list):

    input_flirted_path = InputPath.split('.')[0]+'_flirted.nii.gz'
    input_flirt_path = InputPath.split('.')[0]+'_flirt.mat'
    input_invflirt_path = InputPath.split('.')[0]+'_invflirt.mat'
    input_flirted_synned_path = InputPath.split('.')[0]+'_flirted_synned.nii.gz'
    annotation_invsynned_path = InputPath.split('.')[0]+'_flirted_annotation.nii.gz'
    annotation_invsynned_invflirted_path = InputPath.split('.')[0]+'_annotation.nii.gz'

    # FLIRT input to reference space
    os.system('flirt -in ' + InputPath + ' \
                     -ref ' + template_path + ' \
                     -out ' + input_flirted_path + ' \
                     -omat ' + input_flirt_path + ' \
                     -verbose 1')

    # Load images
    annotation_image = nib.load(annotation_path)
    annotation = annotation_image.get_fdata()
    template_image = nib.load(template_path)
    template = template_image.get_fdata()
    input_flirted_image = nib.load(input_flirted_path)
    input_flirted = input_flirted_image.get_fdata()

    # SyN
    metric = CCMetric(3)
    level_iters = [10, 10, 5, 5, 5]
    print(datetime.datetime.now())
    sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)

    mapping = sdr.optimize(static=template, moving=input_flirted,
                           static_grid2world=template_image.get_qform(), moving_grid2world=input_flirted_image.get_qform())

    input_flirted_synned = mapping.transform(input_flirted)

    annotation_invsynned = mapping.transform_inverse(annotation, interpolation='nearest')

    # save some (inv)synned warped images, mostly to check if everything went well
    input_flirted_synned_image = nib.Nifti1Image(input_flirted_synned, np.eye(4))
    input_flirted_synned_image.set_qform(template_image.get_qform(), code=1)
    input_flirted_synned_image.set_sform(np.eye(4), code=0)
    nib.save(input_flirted_synned_image, input_flirted_synned_path)

    annotation_invsynned_image = nib.Nifti1Image(annotation_invsynned, np.eye(4))
    annotation_invsynned_image.set_qform(input_flirted_image.get_qform(), code=1)
    annotation_invsynned_image.set_sform(np.eye(4), code=0)
    nib.save(annotation_invsynned_image, annotation_invsynned_path)

    # inflirt invsynned annotation to flirted image to get annotation of original image
    os.system('convert_xfm -omat '+input_invflirt_path+' -inverse '+input_flirt_path)
    os.system('flirt -in ' + annotation_invsynned_path + ' \
                     -ref ' + InputPath + ' \
                     -out ' + annotation_invsynned_invflirted_path + ' \
                     -init ' + input_invflirt_path + ' \
                     -applyxfm \
                     -interp nearestneighbour \
                     -verbose 1')
