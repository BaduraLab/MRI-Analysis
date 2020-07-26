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
annotation_path_list = [os.path.join(reference_path, 'atlases', 'Cerebellum', 'Cerebellum-MNIfnirt-prob-1mm.nii.gz'),
                        os.path.join(reference_path, 'atlases', 'Talairach', 'Talairach-labels-1mm.nii.gz')]
template_path = os.path.join(reference_path, 'standard', 'MNI152_T1_1mm_brain.nii.gz')
input_path_list = glob.glob(os.path.join(data_path, '*', '*_reoriented.nii.gz'))



# Loop through inputs
for iInputPath, InputPath in enumerate(input_path_list):
    print(iInputPath)
    print(InputPath)
    print(datetime.datetime.now())

    input_flirted_path = InputPath.split('.')[0]+'_flirted.nii.gz'
    input_flirt_path = InputPath.split('.')[0]+'_flirt.mat'
    input_invflirt_path = InputPath.split('.')[0]+'_invflirt.mat'
    input_flirted_synned_path = InputPath.split('.')[0]+'_flirted_synned.nii.gz'
    annotation_invsynned_invflirted_4D_path = InputPath.split('.')[0]+'_annotation_4D.nii.gz'

    # FLIRT input to reference space
    os.system('flirt -in ' + InputPath + ' \
                     -ref ' + template_path + ' \
                     -out ' + input_flirted_path + ' \
                     -omat ' + input_flirt_path + ' \
                     -verbose 1')

    # Load images
    input_image = nib.load(InputPath)
    template_image = nib.load(template_path)
    template = template_image.get_fdata()
    input_flirted_image = nib.load(input_flirted_path)
    input_flirted = input_flirted_image.get_fdata()

    # SyN
    metric = CCMetric(3)
    level_iters = [10, 10, 5, 5, 5]
    sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)

    mapping = sdr.optimize(static=template, moving=input_flirted,
                           static_grid2world=template_image.get_qform(), moving_grid2world=input_flirted_image.get_qform())

    input_flirted_synned = mapping.transform(input_flirted)

    # save the flirted and synned input, this is the input aligned to the reference elastically
    input_flirted_synned_image = nib.Nifti1Image(input_flirted_synned, np.eye(4))
    input_flirted_synned_image.set_qform(template_image.get_qform(), code=1)
    input_flirted_synned_image.set_sform(np.eye(4), code=0)
    nib.save(input_flirted_synned_image, input_flirted_synned_path)

    # Begin annotation loop
    annotation_invsynned_path_list = list()
    annotation_invsynned_invflirted_path_list = list()
    annotation_invsynned_invflirted_list = list()
    for iAnnotationPath, AnnotationPath in enumerate(annotation_path_list):
        print(iAnnotationPath)

        annotation_4D_image = nib.load(AnnotationPath)
        annotation_4D = annotation_4D_image.get_fdata()
        annotation_4D_shape = annotation_4D.shape
        interpolation_method = 'linear'
        if len(annotation_4D_shape) == 3: # if it is not actually a 4D image, extend it to be 4D
            annotation_4D_shape = annotation_4D_shape + tuple([1])
            annotation_4D = annotation_4D.reshape(annotation_4D_shape)
            interpolation_method = 'nearest'

        # Loop over 4th dimension of annotation, for label annotation this defaults to one iteration (no loop really)
        annotation_invsynned_path_list_perannotation = list()
        annotation_invsynned_invflirted_path_list_perannotation = list()
        annotation_invsynned_invflirted_list_perannotation = list()
        for i4D in range(annotation_4D_shape[3]):
            print(i4D)

            annotation = annotation_4D[:, :, :, i4D]
            annotation_invsynned = mapping.transform_inverse(annotation, interpolation=interpolation_method) # nearest also works for non-probabilistic atlases!

            # Define annotation output paths automatically while saving them to lists
            annotation_name = AnnotationPath.split(os.sep)[-1].split('.')[0]
            annotation_invsynned_path = InputPath.split('.')[0]+'_flirted_annotation_'+annotation_name+'_'+str(i4D)+'.nii.gz'
            annotation_invsynned_invflirted_path = InputPath.split('.')[0]+'_annotation_'+annotation_name+'_'+str(i4D)+'.nii.gz'
            annotation_invsynned_path_list.append(annotation_invsynned_path)
            annotation_invsynned_invflirted_path_list.append(annotation_invsynned_invflirted_path)

            annotation_invsynned_path_list_perannotation.append(annotation_invsynned_path)
            annotation_invsynned_invflirted_path_list_perannotation.append(annotation_invsynned_invflirted_path)

            # save invsynned annotation (this is then the annotation of the flirted image)
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

            # Load annotation of native image and save to list
            annotation_invsynned_invflirted_image = nib.load(annotation_invsynned_invflirted_path)
            annotation_invsynned_invflirted = annotation_invsynned_invflirted_image.get_fdata()
            annotation_invsynned_invflirted_list.append(annotation_invsynned_invflirted)

            annotation_invsynned_invflirted_list_perannotation.append(annotation_invsynned_invflirted)

        # Concatenate in 4D and save with nib
        annotation_invsynned_invflirted_4D_perannotation = np.array(annotation_invsynned_invflirted_list_perannotation).transpose(1, 2, 3, 0)
        annotation_invsynned_invflirted_4D_image_perannotation = nib.Nifti1Image(annotation_invsynned_invflirted_4D_perannotation, np.eye(4))
        annotation_invsynned_invflirted_4D_image_perannotation.set_qform(input_image.get_qform(), code=1)
        annotation_invsynned_invflirted_4D_image_perannotation.set_sform(np.eye(4), code=0)
        annotation_invsynned_invflirted_4D_path_perannotation = InputPath.split('.')[0]+'_annotation_'+annotation_name+'.nii.gz'
        nib.save(annotation_invsynned_invflirted_4D_image_perannotation, annotation_invsynned_invflirted_4D_path_perannotation)
        # End annotation loop, having loaded individual annotations of different atlases

    # Concatenate in 4D and save with nib
    annotation_invsynned_invflirted_4D = np.array(annotation_invsynned_invflirted_list).transpose(1, 2, 3, 0)
    annotation_invsynned_invflirted_4D_image = nib.Nifti1Image(annotation_invsynned_invflirted_4D, np.eye(4))
    annotation_invsynned_invflirted_4D_image.set_qform(input_image.get_qform(), code=1)
    annotation_invsynned_invflirted_4D_image.set_sform(np.eye(4), code=0)
    nib.save(annotation_invsynned_invflirted_4D_image, annotation_invsynned_invflirted_4D_path)

