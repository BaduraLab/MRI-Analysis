#!/usr/bin/python3
""" Linearly and non-linearly transform human subjects to reference space
"""

__author__ = "Enzo Nio"
__version__ = "1.0.0"
__maintainer__ = "Enzo Nio"

import os
import nibabel as nib
import datetime
import numpy as np
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.imwarp import DiffeomorphicMap
from dipy.align.metrics import CCMetric
# from dipy.align.metrics import EMMetric
# from dipy.align.metrics import SSDMetric
import glob
from pathlib import Path
from functions import save_image
from compress_pickle import dump, load



# Define paths
data_path = os.path.join('Data', 'Human', 'Processed')
# processed_path = os.path.join('Data', 'Human', 'Processed_EMMetric')
# processed_path = os.path.join('Data', 'Human', 'Processed_SSDMetric')
Path(data_path).mkdir(exist_ok=True)
reference_path = os.path.join('Data', 'Human', 'Reference')
# annotation_path = os.path.join('atlases', 'Cerebellum', 'Talairach', 'Talairach-labels-1mm.nii.gz')
annotation_path_list = [os.path.join(reference_path,
                                     'atlases',
                                     'Cerebellum',
                                     'Cerebellum-MNIfnirt-prob-1mm_reoriented.nii.gz'),
                        os.path.join(reference_path,
                                     'subcortical',
                                     'prob_atlas_bilateral_reoriented.nii.gz'),
                        os.path.join(reference_path,
                                     'CerebrA',
                                     'mni_icbm152_CerebrA_tal_nlin_sym_09c_reoriented.nii.gz'),
                        os.path.join(reference_path,
                                     'CerebrA',
                                     'mni_icbm152_t1_tal_nlin_sym_09c_mask_reoriented.nii.gz'),
                        os.path.join(reference_path,
                                     'AAN',
                                     'AAN_reoriented.nii.gz'),
                        os.path.join(reference_path,
                                     'allen',
                                     'annotation_full_custom_reoriented.nii.gz')] ###################################### ADD
template_path_list = [os.path.join(reference_path,
                                   'standard',
                                   'MNI152_T1_1mm_brain_reoriented.nii.gz'),
                      os.path.join(reference_path,
                                   'subcortical',
                                   'CIT168_T1w_700um_reoriented.nii.gz'),
                      os.path.join(reference_path,
                                   'CerebrA',
                                   'mni_icbm152_t1_tal_nlin_sym_09c_masked_reoriented.nii.gz'),
                      os.path.join(reference_path,
                                   'CerebrA',
                                   'mni_icbm152_t1_tal_nlin_sym_09c_masked_reoriented.nii.gz'),
                      os.path.join(reference_path,
                                   'standard',
                                   'MNI152_T1_1mm_brain_reoriented.nii.gz'),
                        os.path.join(reference_path,
                                     'allen',
                                     'mni_icbm152_t1_tal_nlin_asym_09b_hires_brain_reoriented.nii.gz')]
annotation_name_list = ['suit',
                        'subcortical',
                        'CerebrA',
                        'mask',
                        'AAN',
                        'allen']
annotation_path_list = [annotation_path_list[-1]]
template_path_list = [template_path_list[-1]]
annotation_name_list = [annotation_name_list[-1]]
# annotation_path_list = [os.path.join(reference_path,
#                                      'subcortical',
#                                      'prob_atlas_bilateral_reoriented.nii.gz')]
# template_path_list = [os.path.join(reference_path,
#                                    'subcortical',
#                                    'CIT168_T1w_700um_reoriented.nii.gz')]
# annotation_name_list = ['subcortical']
input_path_list = glob.glob(os.path.join(data_path, '*', '*_reoriented.nii.gz'))
input_skull_path_list = glob.glob(os.path.join(data_path, '*', '*skull_reoriented.nii.gz'))
input_path_list = list(set(input_path_list) - set(input_skull_path_list))
# input_path_list = [input_path_list[4]]
input_orsuit_path_list = glob.glob(os.path.join(data_path, '*', 'iw_Lobules-SUIT_u_a_*_reoriented_seg1.nii'))



# Define
probability_threshold = [0.2, 0.4, np.nan, np.nan, np.nan]  # Might be changed to list the same length as number of annotations
# probability_threshold = 0.5
def saveImage(image_fdata, path, image_qform_template):
    image = nib.Nifti1Image(image_fdata, image_qform_template.affine, image_qform_template.header)
    image.set_qform(image_qform_template.affine, code=1)
    image.set_sform(image_qform_template.affine, code=0)
    nib.save(image, path)



# Loop through inputs
for iInputPath, InputPath in enumerate(input_path_list):
    print(iInputPath)
    print(InputPath)
    print(datetime.datetime.now())

    input_dirname = os.path.dirname(InputPath)
    input_name = os.path.basename(InputPath).split('_')[0]
    input_noext = os.path.join(input_dirname, input_name)

    # image with skull path
    input_skull_path = os.path.join(input_dirname, input_name+'_skull_reoriented.nii.gz')
    input_skull_flirted_path = os.path.join(input_dirname, input_name+'_skull_reoriented.nii.gz')
    print(input_skull_path)

    annotation_invsynned_invflirted_composite_path = input_noext+'_annotation_composite.nii.gz'

    # Begin annotation loop
    annotation_invsynned_invflirted_image_list_list = list()
    for iAnnotationPath, AnnotationPath in enumerate(annotation_path_list):
        template_path = template_path_list[iAnnotationPath]
        # template_name = template_path.split(os.sep)[-1].split('.')[0]
        template_name = annotation_name_list[iAnnotationPath]
        print(iAnnotationPath)

        input_flirted_path = input_noext+'_flirted_'+template_name+'_'+str(iAnnotationPath)+'.nii.gz'
        input_flirt_path = input_noext+'_flirt_'+template_name+'_'+str(iAnnotationPath)+'.mat'
        input_flirted_syn_path = input_noext+'_flirted_syn_'+template_name+'_'+str(iAnnotationPath)+'.pickle.gz'
        input_invflirt_path = input_noext+'_invflirt_'+template_name+'_'+str(iAnnotationPath)+'.mat'
        input_flirted_synned_path = input_noext+'_flirted_synned_'+template_name+'_'+str(iAnnotationPath)+'.nii.gz'
        annotation_name = annotation_name_list[iAnnotationPath]
        # annotation_name = AnnotationPath.split(os.sep)[-1].split('.')[0].split('_')[0]
        annotation_outdir = Path(input_dirname, annotation_name)
        annotation_outdir.mkdir(exist_ok=True)

        # FLIRT input to reference space
        os.system('flirt -in ' + InputPath + ' \
                         -ref ' + template_path + ' \
                         -out ' + input_flirted_path + ' \
                         -omat ' + input_flirt_path + ' \
                         -verbose 0')

        # Load images
        input_image = nib.load(InputPath)
        template_image = nib.load(template_path)
        template = template_image.get_fdata()
        input_flirted_image = nib.load(input_flirted_path)
        input_flirted = input_flirted_image.get_fdata()
        # if os.path.isfile(input_skull_flirted_path):
        if False: # skull processing off
            print('with skull processing')
            os.system('flirt -in ' + input_skull_path + ' \
                             -ref ' + template_path + ' \
                             -out ' + input_skull_flirted_path + ' \
                             -init ' + input_flirt_path + ' \
                             -verbose 0')
            input_skull_flirted_image = nib.load(input_skull_flirted_path)
            input_skull_flirted = input_skull_flirted_image.get_fdata()
        else:
            print('without skull processing')
            input_skull_flirted = input_flirted



        # SyN
        # metric = SSDMetric(3)
        # metric = EMMetric(3)
        metric = CCMetric(3)
        level_iters = [10, 10, 5, 5, 5]
        sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)

        mapping = sdr.optimize(static=template,
                               moving=input_skull_flirted,
                               static_grid2world=template_image.get_qform(),
                               moving_grid2world=input_flirted_image.get_qform())
        with open(input_flirted_syn_path, 'wb') as f:
            dump([mapping, metric, level_iters, sdr], f, protocol=4, compression='gzip')

        forw_field = mapping.get_forward_field()
        back_field = mapping.get_backward_field()
        forw_SS = np.sum(np.power(forw_field, 2))
        back_SS = np.sum(np.power(back_field, 2))
        dif_SSD = np.sum(np.power(forw_field + back_field, 2))
        dif_SSD_norm = dif_SSD / ((forw_SS + back_SS) / 2)
        print(f'dif_SSD_norm = {dif_SSD_norm}')

        input_flirted_synned = mapping.transform(input_flirted)
        input_flirted_synned_invsynned = mapping.transform_inverse(input_flirted_synned)

        # save the flirted and synned input, this is the input aligned to the reference elastically
        input_flirted_synned_image = nib.Nifti1Image(input_flirted_synned, np.eye(4))
        input_flirted_synned_image.set_qform(template_image.get_qform(), code=1)
        input_flirted_synned_image.set_sform(np.eye(4), code=0)
        nib.save(input_flirted_synned_image, input_flirted_synned_path)
        annotation_invsynned_invflirted_4D_path = input_noext+'_annotation_'+annotation_name+'_4D.nii.gz'
        annotation_invsynned_invflirted_maxprob_path = input_noext+'_annotation_'+annotation_name+'_maxprob.nii.gz'
        annotation_invsynned_invflirted_thresholded_path = input_noext+'_annotation_'+annotation_name+'_thrarg.nii.gz'



        annotation_4D_image = nib.load(AnnotationPath)
        annotation_4D = annotation_4D_image.get_fdata()
        annotation_4D_shape = annotation_4D.shape
        interpolation_method = 'linear'
        interpolation_method_flirt = 'trilinear'
        is3D = len(annotation_4D_shape) == 3
        if is3D:  # if it is not actually a 4D image, extend it to be 4D
            annotation_4D_shape = annotation_4D_shape + tuple([1])
            annotation_4D = annotation_4D.reshape(annotation_4D_shape)
            interpolation_method = 'nearest'
            interpolation_method_flirt = 'nearestneighbour'

        # Loop over 4th dimension of annotation, for label annotation this defaults to one iteration (no loop really)
        annotation_invsynned_invflirted_image_list = list()
        for i4D in range(annotation_4D_shape[3]):
            print('i4D='+str(i4D))

            annotation = annotation_4D[:, :, :, i4D]

            annotation_invsynned = mapping.transform_inverse(annotation, interpolation=interpolation_method) # nearest also works for non-probabilistic atlases!

            # Define annotation output paths automatically while saving them to lists
            annotation_invsynned_path = os.path.join(annotation_outdir,
                                                     input_name+'_flirted_annotation_'+annotation_name+'_'+str(i4D+1)+'.nii.gz')
            annotation_invsynned_invflirted_path = os.path.join(annotation_outdir,
                                                                input_name+'_annotation_'+annotation_name+'_'+str(i4D+1)+'.nii.gz')
            # annotation_invsynned_path_list.append(annotation_invsynned_path)
            # annotation_invsynned_invflirted_path_list.append(annotation_invsynned_invflirted_path)

            # save invsynned annotation (this is then the annotation of the flirted image)
            annotation_invsynned_image = nib.Nifti1Image(annotation_invsynned, np.eye(4))
            annotation_invsynned_image.set_qform(annotation_4D_image.get_qform(), code=1)
            annotation_invsynned_image.set_sform(np.eye(4), code=0)
            nib.save(annotation_invsynned_image, annotation_invsynned_path)

            # invflirt invsynned annotation to flirted image to get annotation of original image
            os.system('convert_xfm -omat '+input_invflirt_path+' -inverse '+input_flirt_path)
            os.system('flirt -in ' + annotation_invsynned_path + ' \
                             -ref ' + InputPath + ' \
                             -out ' + annotation_invsynned_invflirted_path + ' \
                             -init ' + input_invflirt_path + ' \
                             -applyxfm \
                             -interp ' + interpolation_method_flirt + ' \
                             -verbose 0')

            # Load annotation of native image and save to list
            annotation_invsynned_invflirted_image = nib.load(annotation_invsynned_invflirted_path)
            annotation_invsynned_invflirted_image_list.append(annotation_invsynned_invflirted_image)

        # Concatenate in 4D and save with nib
        annotation_invsynned_invflirted_4D_image = nib.concat_images(annotation_invsynned_invflirted_image_list)
        nib.save(annotation_invsynned_invflirted_4D_image, annotation_invsynned_invflirted_4D_path)

        # If originally 4D, create maximum probability image
        if not is3D:
            annotation_invsynned_invflirted_4D = annotation_invsynned_invflirted_4D_image.get_fdata()
            annotation_invsynned_invflirted_maxprob = np.max(annotation_invsynned_invflirted_4D, axis=3)
            annotation_invsynned_invflirted_maxprob = annotation_invsynned_invflirted_maxprob / \
                                                      np.max(annotation_invsynned_invflirted_maxprob)
            annotation_invsynned_invflirted_argprob = np.argmax(annotation_invsynned_invflirted_4D, axis=3)+1
            annotation_invsynned_invflirted_thrprob = annotation_invsynned_invflirted_maxprob > probability_threshold[iAnnotationPath]
            annotation_invsynned_invflirted = annotation_invsynned_invflirted_argprob \
                                            * annotation_invsynned_invflirted_thrprob

            saveImage(image_fdata=annotation_invsynned_invflirted_maxprob,
                      path=annotation_invsynned_invflirted_maxprob_path,
                      image_qform_template=annotation_invsynned_invflirted_image)
            saveImage(image_fdata=annotation_invsynned_invflirted_thrprob.astype(int),
                      path=annotation_invsynned_invflirted_maxprob_path.split('.')[0]+'_thresholded.nii.gz',
                      image_qform_template=annotation_invsynned_invflirted_image)
            saveImage(image_fdata=annotation_invsynned_invflirted.astype(np.int16),
                      path=annotation_invsynned_invflirted_thresholded_path,
                      image_qform_template=annotation_invsynned_invflirted_image)
        else:
            annotation_invsynned_invflirted = annotation_invsynned_invflirted_4D_image.get_fdata()
            saveImage(image_fdata=annotation_invsynned_invflirted.astype(np.int16),
                      path=annotation_invsynned_invflirted_thresholded_path,
                      image_qform_template=annotation_invsynned_invflirted_4D_image)


        # Add list of slices of 4D images to list for creation of composite image
        annotation_invsynned_invflirted_image_list_list.append(annotation_invsynned_invflirted_image_list)
        # End annotation loop

    # Concatenate in 4D (different atlases) and save with nib
    annotation_invsynned_invflirted_composite_image = \
        nib.concat_images(sum(annotation_invsynned_invflirted_image_list_list, []))
    nib.save(annotation_invsynned_invflirted_composite_image,
             annotation_invsynned_invflirted_composite_path)



# Copy images processed by MATLAB SUIT to files with your naming scheme
for Path in input_orsuit_path_list:
    orsuit_image = nib.load(Path)

    subject = Path.split(os.sep)[-2]

    # Define output paths
    orsuit_adjusted_path = os.path.join(data_path, subject, subject+'_annotation_orsuit_thrarg.nii.gz')
    orsuit_adjusted_gray_path = os.path.join(data_path, subject, subject+'_annotation_orsuit_thrarg_gray.nii.gz')
    # orsuit_adjusted_1_path = os.path.join(data_path, subject, subject+'_annotation_orsuit_thrarg_adjusted_1.nii.gz')
    # orsuit_adjusted_2_path = os.path.join(data_path, subject, subject+'_annotation_orsuit_thrarg_adjusted_2.nii.gz')

    nib.save(orsuit_image, orsuit_adjusted_path)
    print(orsuit_adjusted_path)
    orsuit = orsuit_image.get_fdata()
    orsuit_gray = np.isin(orsuit, np.arange(28)+1)
    saveImage(orsuit_gray, orsuit_adjusted_gray_path, orsuit_image)
    print(orsuit_adjusted_gray_path)
    # if not os.path.isfile(orsuit_adjusted_1_path):
    #     nib.save(orsuit_image, orsuit_adjusted_1_path)
    #     print(orsuit_adjusted_1_path)
    # if not os.path.isfile(orsuit_adjusted_2_path):
    #     nib.save(orsuit_image, orsuit_adjusted_2_path)
    #     print(orsuit_adjusted_2_path)