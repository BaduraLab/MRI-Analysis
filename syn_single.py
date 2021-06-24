# Import
import numpy as np
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.imwarp import DiffeomorphicMap
from dipy.align.metrics import CCMetric
from dipy.core.gradients import gradient_table
from dipy.data import get_fnames
from dipy.io.image import load_nifti, save_nifti
from dipy.io.gradients import read_bvals_bvecs
import os.path
from dipy.viz import regtools
import scipy.io
import glob
import nibabel as nib
from scipy import spatial
import datetime
# import plotly.graph_objects as go

# Define
allen_resolution = 25
# reference_path = '/usr/local/fsl/data/standard/allen_new'
reference_path = os.path.join('Data', 'Mouse', 'Reference')
data_path = os.path.join('Data', 'Mouse', 'Processed')
reference_image_path = os.path.join(reference_path, 'average_template_'+str(allen_resolution)+'_to_AMBMC_flirted.nii.gz')
annotation_image_path = os.path.join(reference_path, 'annotation_'+str(allen_resolution)+'_to_AMBMC_flirted.nii.gz')

# data_path = '/home/enzo/Desktop/Data/Mouse/Processed_New'
mouse_path_list = glob.glob(os.path.join(data_path, '*'))
# mouse_path_list = [s.split('/')[-2] for s in mouse_path_list]
# subjects = subjects[2::] # remove WT_50, which was already synned with a resolution of 25 micrometers

# for iSubject in range(len(subjects)):
mouse_path_list = [mouse_path_list[-3]]
for iMousePath, MousePath in enumerate(mouse_path_list):
    start_time = datetime.datetime.now()
    print(iMousePath)
    print(iMousePath)
    print(MousePath)
    print(f'Start time = {start_time}')

    # Define per suject
    mouse_string = MousePath.split(os.sep)[-1]
    # subject = subjects[iSubject]
    # input_path = os.path.join(data_path, mouse_string)
    FLIRT_folder_path = os.path.join(MousePath, 'FLIRT')
    # mask_path = glob.glob(FLIRT_folder_path + '/*.mask.nii.gz')[0]
    # FLIRT_path = glob.glob(FLIRT_folder_path + '/*warpaffine.mat')[0]
    # FLIRT_inverted_path = glob.glob(FLIRT_folder_path + '/*warpaffine_inverted.mat')[0]
    input_image_path = os.path.join(MousePath, mouse_string+'_reoriented.nii.gz')
    mask_path = os.path.join(MousePath, mouse_string + '_mask_t=500_v=380_k=6.mask.nii.gz')
    mouse_masked_path = os.path.join(MousePath, mouse_string+'_masked.nii.gz')
    mouse_masked_flirted_path = os.path.join(MousePath, mouse_string+'_flirted.nii.gz')
    mouse_masked_flirt_path = os.path.join(MousePath, mouse_string+'_flirt.mat')
    mouse_masked_invflirt_path = os.path.join(MousePath, mouse_string+'_invflirt.mat')
    input_original_image_path = os.path.join(MousePath, mouse_string+'.nii.gz')
    input_invwarped_path = os.path.join(MousePath, mouse_string+'_synned.nii.gz')
    annotation_invwarped_path = os.path.join(MousePath, 'allen_annotation_invsynned_to_'+mouse_string+'_flirted.nii.gz')
    annotation_subject_path = os.path.join(MousePath, 'allen_annotation_invsynned_to_'+mouse_string+'.nii.gz')
    annotation_subject_adjusted_path = os.path.join(MousePath, 'allen_annotation_invsynned_to_'+mouse_string+'_adjusted.nii.gz')

    # Load images
    reference_image = nib.load(reference_image_path)
    reference = reference_image.get_fdata()
    annotation_image = nib.load(annotation_image_path)
    annotation = annotation_image.get_fdata()
    input_image = nib.load(input_image_path)
    input = input_image.get_fdata()



    # Mask mouse image
    mask_image = nib.load(mask_path)
    mask = mask_image.get_fdata()
    mask = mask / np.max(mask)
    print('mask image')
    mouse_masked = input * mask
    mouse_masked_image = nib.Nifti1Image(mouse_masked, input_image.affine, input_image.header)
    nib.save(mouse_masked_image, mouse_masked_path)

    # # FLIRT subject to reference
    # print('FLIRT rigid start')
    # os.system('flirt -in ' + mouse_masked_path + ' \
    #                  -ref ' + reference_template_path + ' \
    #                  -out ' + mouse_masked_flirtedRigid_path + ' \
    #                  -omat ' + mouse_masked_flirtRigid_path + ' \
    #                  -dof ' + '6' + ' \
    #                  -verbose 0')    # FLIRT subject to reference

    # FLIRT rigidly transformed subject to reference
    print('FLIRT affine start')
    os.system('flirt -in ' + mouse_masked_path + ' \
                     -ref ' + reference_image_path + ' \
                     -out ' + mouse_masked_flirted_path + ' \
                     -omat ' + mouse_masked_flirt_path + ' \
                     -verbose 0')
    mouse_masked_flirted_image = nib.load(mouse_masked_flirted_path)
    mouse_masked_flirted = mouse_masked_flirted_image.get_fdata()





    # syn
    metric = CCMetric(3)
    level_iters = [10, 10, 5, 5, 5]
    print(iMousePath)
    print(datetime.datetime.now())
    sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)

    mapping = sdr.optimize(static=reference, moving=mouse_masked_flirted,
                           static_grid2world=reference_image.get_qform(), moving_grid2world=input_image.get_qform())

    warped_image = mapping.transform(mouse_masked_flirted)

    warped_reference = mapping.transform_inverse(reference)

    warped_annotation = mapping.transform_inverse(annotation, interpolation='nearest')



    # save warps
    input_invwarped_image = nib.Nifti1Image(warped_image, np.eye(4))
    input_invwarped_image.set_qform(reference_image.get_qform(), code=1)
    input_invwarped_image.set_sform(np.eye(4), code=0)
    nib.save(input_invwarped_image, input_invwarped_path)

    annotation_invwarped_image = nib.Nifti1Image(warped_annotation, np.eye(4))
    annotation_invwarped_image.set_qform(reference_image.get_qform(), code=1)
    annotation_invwarped_image.set_sform(np.eye(4), code=0)
    nib.save(annotation_invwarped_image, annotation_invwarped_path)



    # # evaluate warps
    # X, Y, Z = np.mgrid[0:input.shape[0]:1, 650:651:1, 0:input.shape[2]:1]
    #
    # fig1 = go.Figure(data=go.Volume(
    #     x=X.flatten(),
    #     y=Y.flatten(),
    #     z=Z.flatten(),
    #     value=input[0:input.shape[0]:1, 650:651:1, 0:input.shape[2]:1].flatten(),
    #     opacity=1,  # needs to be small to see through all surfaces
    #     surface_count=17,  # needs to be a large number for good volume rendering
    # ))
    # fig1.show()
    #
    # X, Y, Z = np.mgrid[0:reference.shape[0]:1, 650:651:1, 0:reference.shape[2]:1]
    #
    # fig2 = go.Figure(data=go.Volume(
    #     x=X.flatten(),
    #     y=Y.flatten(),
    #     z=Z.flatten(),
    #     value=reference[0:reference.shape[0]:1, 650:651:1, 0:reference.shape[2]:1].flatten(),
    #     opacity=1,  # needs to be small to see through all surfaces
    #     surface_count=17,  # needs to be a large number for good volume rendering
    # ))
    # fig2.show()
    #
    # X, Y, Z = np.mgrid[0:warped_image.shape[0]:1, 650:651:1, 0:warped_image.shape[2]:1]
    #
    # fig3 = go.Figure(data=go.Volume(
    #     x=X.flatten(),
    #     y=Y.flatten(),
    #     z=Z.flatten(),
    #     value=warped_image[0:warped_image.shape[0]:1, 650:651:1, 0:warped_image.shape[2]:1].flatten(),
    #     opacity=1,  # needs to be small to see through all surfaces
    #     surface_count=17,  # needs to be a large number for good volume rendering
    # ))
    # fig3.show()
    #
    # X, Y, Z = np.mgrid[0:warped_annotation.shape[0]:1, 650:651:1, 0:warped_annotation.shape[2]:1]
    #
    # fig4 = go.Figure(data=go.Volume(
    #     x=X.flatten(),
    #     y=Y.flatten(),
    #     z=Z.flatten(),
    #     value=warped_annotation[0:warped_annotation.shape[0]:1, 650:651:1, 0:warped_annotation.shape[2]:1].flatten(),
    #     opacity=1,  # needs to be small to see through all surfaces
    #     surface_count=17,  # needs to be a large number for good volume rendering
    # ))
    # fig4.show()

    # Invert FLIRT warped annotation back to subject space
    print('inverse affine FLIRT')
    os.system('convert_xfm -omat '+mouse_masked_invflirt_path+' -inverse '+mouse_masked_flirt_path)
    os.system('flirt -in ' + annotation_invwarped_path + ' \
                     -ref ' + input_original_image_path + ' \
                     -out ' + annotation_subject_path + ' \
                     -init ' + mouse_masked_invflirt_path + ' \
                     -applyxfm \
                     -interp nearestneighbour \
                     -verbose 1')



    ## Adjust synned annotation
    # load images
    annotation_subject_image = nib.load(annotation_subject_path)
    annotation_subject = annotation_subject_image.get_fdata()
    input_original_image = nib.load(input_original_image_path)
    input_original = input_original_image.get_fdata()

    # make annotation 0 everywhere outside of mask
    annotation_subject_adjusted = annotation_subject * mask

    # save adjusted invsynned
    annotation_subject_adjusted_image = nib.Nifti1Image(annotation_subject_adjusted, np.eye(4))
    annotation_subject_adjusted_image.set_qform(input_original_image.get_qform(), code=1)
    annotation_subject_adjusted_image.set_sform(np.eye(4), code=0)
    nib.save(annotation_subject_adjusted_image, annotation_subject_adjusted_path)

    # X, Y, Z = np.mgrid[0:annotation_subject.shape[0]:1, 0:annotation_subject.shape[1]:1, 0:annotation_subject.shape[2]:1]
    #
    # # points which are unannotated and inside mask
    # unannotated = annotation_subject == 0
    # unannotated_inside = np.logical_and(unannotated, mask != 0)
    # unannotated_inside_points = np.vstack((X[unannotated_inside], Y[unannotated_inside], Z[unannotated_inside])).transpose()
    #
    # # points which have annotation (also outside mask)
    # annotated = np.logical_not(unannotated)
    # annotated_points = np.vstack((X[annotated], Y[annotated], Z[annotated])).transpose()
    #
    # # create KDtree
    # tree = spatial.KDTree(annotated_points)
    #
    # # for each point without annotation inside mask, find nearest point with annotation within or outside mask
    # for iUnannotatedInside in range(len(unannotated_inside_points)):
    #     unannotated_inside_index = tuple(unannotated_inside_points[iUnannotatedInside])
    #     closest_annotated_index = tuple(annotated_points[tree.query(unannotated_inside_index)[1]])
    #     annotation_subject_adjusted[unannotated_inside_index] = annotation_subject[closest_annotated_index]
    #
    # # save adjusted invsynned
    # annotation_subject_adjusted_image = nib.Nifti1Image(annotation_subject_adjusted, np.eye(4))
    # annotation_subject_adjusted_image.set_qform(input_original_image.get_qform(), code=1)
    # annotation_subject_adjusted_image.set_sform(np.eye(4), code=0)
    # nib.save(annotation_subject_adjusted_image, annotation_subject_adjusted_path)
