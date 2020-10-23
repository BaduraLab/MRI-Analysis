import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
# import fsl.wrappers
import csv
from scipy.stats import ttest_ind
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.imwarp import DiffeomorphicMap
from dipy.align.metrics import CCMetric
import datetime
# import SimpleITK as sitk
from compress_pickle import dump, load

# Define
data_path = os.path.join('Data', 'Mouse', 'Processed')
mouse_path_list = glob.glob(os.path.join(data_path, '*'))
reference_path = os.path.join('Data', 'Mouse', 'Reference')
# average_template_50_to_AMBMC_flirted.nii.gz
# reference_template_path = os.path.join(reference_path, 'average_template_50_reoriented.nii.gz')
# reference_annotation_path = os.path.join(reference_path, 'annotation_50_reoriented.nii.gz')
reference_template_path = os.path.join(reference_path, 'average_template_25_reoriented_flirted_cropped.nii.gz')
reference_annotation_path = os.path.join(reference_path, 'annotation_25_reoriented_flirted_cropped.nii.gz')

# Loop through mice
for iMousePath, MousePath in enumerate(mouse_path_list):
    print(iMousePath)
    print(MousePath)
    print(datetime.datetime.now())

    # Define mouse paths
    mouse_string = MousePath.split(os.sep)[-1]
    mouse_path = os.path.join(MousePath, mouse_string + '_reoriented.nii.gz')
    mask_path = os.path.join(MousePath, mouse_string + '_mask_t=500_v=380_k=6.mask.nii.gz')
    mouse_masked_path = os.path.join(MousePath, mouse_string + '_masked.nii.gz')
    mouse_masked_translated_path = os.path.join(MousePath, mouse_string + '_translated.nii.gz')
    mouse_masked_flirtedRigid_path = os.path.join(MousePath, mouse_string + '_flirtedRigid.nii.gz')
    mouse_masked_flirted_path = os.path.join(MousePath, mouse_string + '_flirted.nii.gz')
    mouse_masked_flirtRigid_path = os.path.join(MousePath, mouse_string + '_flirtRigid.mat')
    mouse_masked_flirt_path = os.path.join(MousePath, mouse_string + '_flirt.mat')
    mouse_masked_flirted_syn_path = os.path.join(MousePath, mouse_string + '_flirted_syn.pickle.gz')
    mouse_masked_invflirt_path = os.path.join(MousePath, mouse_string + '_invflirt.mat')
    mouse_masked_invflirtRigid_path = os.path.join(MousePath, mouse_string + '_invflirtRigid.mat')
    mouse_masked_flirted_synned_path = os.path.join(MousePath, mouse_string + '_flirted_synned.nii.gz')
    reference_annotation_invsynned_path = os.path.join(MousePath, mouse_string + '_annotation_flirted.nii.gz')
    reference_annotation_invsynned_invflirted_path = os.path.join(MousePath, mouse_string + '_annotation_flirtedRigid.nii.gz')
    reference_annotation_invsynned_invflirted_invflirtedRigid_path = os.path.join(MousePath, mouse_string + '_annotation.nii.gz')
    forw_field_path = os.path.join(MousePath, mouse_string + '_annotation.nii.gz')
    back_field_path = os.path.join(MousePath, mouse_string + '_annotation.nii.gz')

    # Load images
    print('loading images')
    mouse_image = nib.load(mouse_path)
    mouse = mouse_image.get_fdata()
    mask_image = nib.load(mask_path)
    mask = mask_image.get_fdata()
    reference_template_image = nib.load(reference_template_path)
    reference_template = reference_template_image.get_fdata()
    reference_annotation_image = nib.load(reference_annotation_path)
    reference_annotation = reference_annotation_image.get_fdata()

    # Mask mouse image
    print('mask image')
    mask = mask / np.max(mask)
    mouse_masked = mouse * mask
    mouse_masked_image = nib.Nifti1Image(mouse_masked, mouse_image.affine, mouse_image.header)
    nib.save(mouse_masked_image, mouse_masked_path)

    # FLIRT subject to reference
    print('FLIRT rigid start')
    os.system('flirt -in ' + mouse_masked_path + ' \
                     -ref ' + reference_template_path + ' \
                     -out ' + mouse_masked_flirtedRigid_path + ' \
                     -omat ' + mouse_masked_flirtRigid_path + ' \
                     -dof ' + '6' + ' \
                     -verbose 0')    # FLIRT subject to reference

    # FLIRT rigidly transformed subject to reference
    print('FLIRT affine start')
    os.system('flirt -in ' + mouse_masked_flirtedRigid_path + ' \
                     -ref ' + reference_template_path + ' \
                     -out ' + mouse_masked_flirted_path + ' \
                     -omat ' + mouse_masked_flirt_path + ' \
                     -verbose 0')
    mouse_masked_flirted_image = nib.load(mouse_masked_flirted_path)
    mouse_masked_flirted = mouse_masked_flirted_image.get_fdata()

    # SyN flirted images to reference
    print('SyN')
    metric = CCMetric(3)
    level_iters = [10, 10, 5, 5, 5]
    sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)
    mapping = sdr.optimize(static=reference_template,
                           moving=mouse_masked_flirted,
                           static_grid2world=reference_template_image.get_qform(),
                           moving_grid2world=mouse_masked_flirted_image.get_qform())
    with open(mouse_masked_flirted_syn_path, 'wb') as f:
        dump([mapping, metric, level_iters, sdr], f, protocol=4, compression='gzip')
    # with open(mouse_masked_syn_path, 'rb') as f:
    #     [mapping, metric, level_iters, sdr] = load(f)

    forw_field = mapping.get_forward_field()
    back_field = mapping.get_backward_field()
    forw_SS = np.sum(np.power(forw_field, 2))
    back_SS = np.sum(np.power(back_field, 2))
    dif_SSD = np.sum(np.power(forw_field + back_field, 2))
    dif_SSD_norm = dif_SSD / ((forw_SS + back_SS) / 2)
    print(f'dif_SSD_norm = {dif_SSD_norm}')

    # forw_field_image = nib.Nifti1Image(forw_field,
    #                                    mouse_masked_flirted_image.affine)
    # nib.save(forw_field_image, forw_field_path)
    # back_field_image = nib.Nifti1Image(back_field,
    #                                    mouse_masked_flirted_image.affine)
    # nib.save(back_field_image, back_field_path)


    mouse_masked_flirted_synned = mapping.transform(mouse_masked_flirted)
    mouse_masked_flirted_synned_image = nib.Nifti1Image(mouse_masked_flirted_synned,
                                                                 mouse_masked_flirted_image.affine,
                                                        mouse_masked_flirted_image.header)
    nib.save(mouse_masked_flirted_synned_image, mouse_masked_flirted_synned_path)

    # Inverse SyN reference to subject space
    print('inverse SyN')
    reference_annotation_invsynned = mapping.transform_inverse(reference_annotation, interpolation='nearest')
    reference_annotation_invsynned_image = nib.Nifti1Image(reference_annotation_invsynned,
                                                           mouse_masked_flirted_image.affine,
                                                           mouse_masked_flirted_image.header)
    nib.save(reference_annotation_invsynned_image, reference_annotation_invsynned_path)

    # invflirt invsynned annotation to flirted image to get annotation of original image
    print('inverse affine FLIRT')
    os.system('convert_xfm -omat '+mouse_masked_invflirt_path+' -inverse '+mouse_masked_flirt_path)
    os.system('flirt -in ' + reference_annotation_invsynned_path + ' \
                     -ref ' + mouse_path + ' \
                     -out ' + reference_annotation_invsynned_invflirted_path + ' \
                     -init ' + mouse_masked_invflirt_path + ' \
                     -applyxfm \
                     -interp nearestneighbour \
                     -verbose 1')

    # invflirt invsynned annotation to flirted image to get annotation of original image
    print('inverse rigid FLIRT')
    os.system('convert_xfm -omat '+mouse_masked_invflirtRigid_path+' -inverse '+mouse_masked_flirtRigid_path)
    os.system('flirt -in ' + reference_annotation_invsynned_invflirted_path + ' \
                     -ref ' + mouse_path + ' \
                     -out ' + reference_annotation_invsynned_invflirted_invflirtedRigid_path + ' \
                     -init ' + mouse_masked_invflirtRigid_path + ' \
                     -applyxfm \
                     -interp nearestneighbour \
                     -verbose 1')



    # invflirt annotation!
    # reference_annotation_invsynned_invflirted_path

    # # FLIRT images to reference
    # elastixImageFilter = sitk.ElastixImageFilter()
    # elastixImageFilter.SetFixedImage(sitk.ReadImage(reference_template_path))
    # elastixImageFilter.SetMovingImage(sitk.ReadImage(mouse_path))
    # elastixImageFilter.SetParameterMap(sitk.GetDefaultParameterMap('translation'))
    # elastixImageFilter.Execute()
    # mouse_masked_translated = elastixImageFilter.GetResultImage()
    # sitk.WriteImage(mouse_masked_translated, mouse_masked_translated_path)
    #
    # elastixImageFilter.SetMovingImage(sitk.ReadImage(mouse_masked_translated_path))
    # elastixImageFilter.SetParameterMap(sitk.GetDefaultParameterMap('affine'))
    # elastixImageFilter.Execute()
    # mouse_masked_flirted = elastixImageFilter.GetResultImage()
    # sitk.WriteImage(mouse_masked_flirted, mouse_masked_flirted_path)
    #
    # mouse_masked_flirted_image = nib.load(mouse_masked_flirted_path)
    # mouse_masked_flirted = mouse_image.get_fdata()

    # fsl.wrappers.flirt(src=mouse_masked_path,
    #                    ref=reference_template_path,
    #                    omat=mouse_masked_flirt_path,
    #                    out=mouse_masked_flirted_path,
    #                    verbose=1)
    # mouse_masked_flirted_image = nib.load(mouse_masked_flirted_path)
    # mouse_masked_flirted = mouse_masked_flirted_image.get_fdata()
    # flirt - in $image_inmasked \
    #             - ref $image_reference \
    #                    - omat $image_warpaffine \
    #                            - out $image_inmasked_flirted \
    #                                   - verbose
    # 1
    # flirt - in $image \
    #             - ref $image_reference \
    #                    - out $image_flirted \
    #                           - init $image_warpaffine \
    #                                   - applyxfm \
    #                                   - verbose 1