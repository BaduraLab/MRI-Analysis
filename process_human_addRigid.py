import os
import numpy as np
import glob
import datetime
from functions import zeroPadImage
import nibabel as nib
from functions import save_image

### Change up the paths

# Define
data_path = os.path.join('Data', 'Human', 'Processed')
subject_path_list = glob.glob(os.path.join(data_path, '*'))
reference_path = os.path.join('Data', 'Human', 'Reference')
reference_template_path = os.path.join(reference_path, 'standard', 'MNI152_T1_0.5mm.nii.gz')
# reference_template_image = nib.load(reference_template_path)
# reference_template = reference_template_image.get_fdata()
# average_template_50_to_AMBMC_flirted.nii.gz
# reference_template_path = os.path.join(reference_path, 'average_template_50_reoriented.nii.gz')
# reference_annotation_path = os.path.join(reference_path, 'annotation_50_reoriented.nii.gz')
# reference_template_path = os.path.join(reference_path, 'average_template_25_reoriented_flirted.nii.gz')
# reference_template_path = os.path.join(reference_path, 'annotation_25_to_AMBMC_flirted.nii.gz')
# reference_annotation_path = os.path.join(reference_path, 'annotation_25_reoriented_flirted.nii.gz')

# Loop through mice
# subject_path_list = subject_path_list[0:1]
for iSubjectPath, SubjectPath in enumerate(subject_path_list):
    print(iSubjectPath)
    print(iSubjectPath)
    print(SubjectPath)
    print(datetime.datetime.now())



    # Define subject paths
    subject_string = SubjectPath.split(os.sep)[-1]
    template_path = os.path.join(SubjectPath, subject_string + '_reoriented.nii.gz')

    # reference_template_path = os.path.join(SubjectPath, subject_string + '_flirted_mask_3.nii.gz')
    # reference_template_image = nib.load(reference_template_path)
    # reference_template = reference_template_image.get_fdata()

    # reference_template_zeropadded, cropped_index = zeroPadImage(reference_template, reference_template, 0.75)
    # reference_template_zeropadded_path = os.path.join(SubjectPath, subject_string + '_flirted_mask_3_zeropadded.nii.gz')
    # save_image(reference_template_zeropadded, reference_template_image, reference_template_zeropadded_path)

    # template_flirtedAffine_path = os.path.join(SubjectPath, subject_string + '_flirtedAffine.nii.gz')
    # flirtAffine_path = os.path.join(SubjectPath, 'flirtAffine.mat')

    annotation_path_list = [
        os.path.join(SubjectPath, subject_string + '_annotation_orsuit_thrarg_adjusted_lobular.nii.gz'),
        os.path.join(SubjectPath, subject_string + '_annotation_subcortical_thrarg.nii.gz'),
        os.path.join(SubjectPath, subject_string + '_annotation_CerebrA_thrarg_adjusted.nii.gz'),
        os.path.join(SubjectPath, subject_string + '_annotation_mask_thrarg.nii.gz'),
        os.path.join(SubjectPath, subject_string + '_annotation_AAN_thrarg.nii.gz')]

    template_flirtedRigid_path = os.path.join(SubjectPath, subject_string + '_flirtedRigid.nii.gz')
    # annotation_flirtedRigid_path = os.path.join(SubjectPath, 'allen_annotation_invsynned_to_' + subject_string + '_adjusted_cerebellum_lobular_flirtedRigid.nii.gz')
    flirtRigid_path = os.path.join(SubjectPath, 'flirtRigid.mat')



    # Affine FLIRT inmasked subject to reference
    print('FLIRT rigid start')
    os.system('flirt -in ' + template_path + ' \
                     -ref ' + reference_template_path + ' \
                     -out ' + template_flirtedRigid_path + ' \
                     -omat ' + flirtRigid_path + ' \
                     -dof ' + '6' + ' \
                     -verbose 0')    # FLIRT subject to reference
    #                  # -dof ' + '6' + ' \

    # # Load affine transformation and separate translation and rotation parameters to
    # flirtAffine = np.loadtxt(flirtAffine_path)
    # trans_matrix = np.identity(4)
    # trans_matrix[0, 3] = flirtAffine[0, 3]
    # trans_matrix[1, 3] = flirtAffine[1, 3]
    # trans_matrix[2, 3] = flirtAffine[2, 3] # translation matrix defined
    # [rotat_matrix, R] = np.linalg.qr(flirtAffine[0:3,0:3]) # rotation matrix defined
    # rigid_matrix = np.identity(4)
    # rigid_matrix[0:3, 0:3] = rotat_matrix
    # rigid_matrix = np.matmul(trans_matrix, rigid_matrix)
    # np.savetxt(flirtRigid_path, rigid_matrix, delimiter='  ', fmt='%.10f')


    # reference_template_path = os.path.join(reference_path,
    #                                        'subcortical',
    #                                        'CIT168_T1w_700um_reoriented.nii.gz')

    # Apply rigid FLIRT to template
    print('FLIRT rigid apply to template start')
    os.system('flirt -in ' + template_path + ' \
                     -ref ' + reference_template_path + ' \
                     -out ' + template_flirtedRigid_path + ' \
                     -init ' + flirtRigid_path + ' \
                     -applyxfm \
                     -interp trilinear \
                     -verbose 1')
    template_flirtedRigid_image = nib.load(template_flirtedRigid_path)
    template_flirtedRigid = template_flirtedRigid_image.get_fdata()
    template_flirtedRigid_zeropadded, cropped_index = zeroPadImage(template_flirtedRigid, template_flirtedRigid, 0.1)
    template_flirtedRigid_zeropadded_path = os.path.join(SubjectPath, subject_string + '_flirtedRigid_cropped.nii.gz')
    save_image(template_flirtedRigid_zeropadded, template_flirtedRigid_image, template_flirtedRigid_zeropadded_path)

    # Apply rigid FLIRT to annotation
    for annotation_path in annotation_path_list:
        annotation_flirtedRigid_path = annotation_path.split('.')[0] + '_flirtedRigid.nii.gz'
        annotation_flirtedRigid_zeropadded_path = annotation_path.split('.')[0] + '_flirtedRigid_cropped.nii.gz'
        print('FLIRT rigid apply to annotation start')
        os.system('flirt -in ' + annotation_path + ' \
                         -ref ' + reference_template_path + ' \
                         -out ' + annotation_flirtedRigid_path + ' \
                         -init ' + flirtRigid_path + ' \
                         -applyxfm \
                         -interp nearestneighbour \
                         -verbose 1')
        annotation_flirtedRigid_image = nib.load(annotation_flirtedRigid_path)
        annotation_flirtedRigid = annotation_flirtedRigid_image.get_fdata()
        annotation_flirtedRigid_zeropadded, cropped_index = zeroPadImage(annotation_flirtedRigid, template_flirtedRigid, 0.1)
        save_image(annotation_flirtedRigid_zeropadded, template_flirtedRigid_image, annotation_flirtedRigid_zeropadded_path)
