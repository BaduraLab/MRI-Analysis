import os
import numpy as np
import glob
import datetime

# Define
data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
mouse_path_list = glob.glob(os.path.join(data_path, '*'))
reference_path = os.path.join('Data', 'Mouse', 'Reference')
# average_template_50_to_AMBMC_flirted.nii.gz
# reference_template_path = os.path.join(reference_path, 'average_template_50_reoriented.nii.gz')
# reference_annotation_path = os.path.join(reference_path, 'annotation_50_reoriented.nii.gz')
reference_template_path = os.path.join(reference_path, 'average_template_25_reoriented_flirted.nii.gz')
# reference_annotation_path = os.path.join(reference_path, 'annotation_25_reoriented_flirted.nii.gz')

# Loop through mice
mouse_path_list = mouse_path_list[0:1]
for iMousePath, MousePath in enumerate(mouse_path_list):
    print(iMousePath)
    print(iMousePath)
    print(MousePath)
    print(datetime.datetime.now())

    # Define mouse paths
    mouse_string = MousePath.split(os.sep)[-1]
    template_path = os.path.join(MousePath, 'FLIRT',  mouse_string + '_inmasked.nii.gz')
    annotation_path = os.path.join(MousePath, 'allen_annotation_invsynned_to_' + mouse_string + '_adjusted_cerebellum_lobular.nii.gz')
    # template_flirtedAffine_path = os.path.join(MousePath, mouse_string + '_inmasked_flirtedAffine.nii.gz')
    flirtAffine_path = os.path.join(MousePath, 'FLIRT', mouse_string + '_to_allen_model_warpaffine.mat')
    template_flirtedRigid_path = os.path.join(MousePath, mouse_string + '_inmasked_flirtedRigid.nii.gz')
    annotation_flirtedRigid_path = os.path.join(MousePath, 'allen_annotation_invsynned_to_' + mouse_string + '_adjusted_cerebellum_lobular_flirtedAffine.nii.gz')
    flirtRigid_path = os.path.join(MousePath, mouse_string + '_flirtRigid.mat')

    # # Affine FLIRT inmasked subject to reference
    # print('FLIRT affine start')
    # os.system('flirt -in ' + template_path + ' \
    #                  -ref ' + reference_template_path + ' \
    #                  -out ' + template_flirtedAffine_path + ' \
    #                  -omat ' + flirtAffine_path + ' \
    #                  -verbose 0')    # FLIRT subject to reference
    #                  # -dof ' + '6' + ' \

    # Load affine transformation and separate translation and rotation parameters to
    flirtAffine = np.loadtxt(flirtAffine_path)
    trans_matrix = np.identity(4)
    trans_matrix[0, 3] = flirtAffine[0, 3]
    trans_matrix[1, 3] = flirtAffine[1, 3]
    trans_matrix[2, 3] = flirtAffine[2, 3] # translation matrix defined
    [rotat_matrix, R] = np.linalg.qr(flirtAffine[0:3,0:3]) # rotation matrix defined
    rigid_matrix = trans_matrix.copy()
    rigid_matrix[0:3, 0:3] = rotat_matrix
    np.savetxt(flirtRigid_path, rigid_matrix, delimiter='  ', fmt='%.10f')

    # Apply rigid FLIRT to template
    print('FLIRT affine apply to annotation start')
    os.system('flirt -in ' + template_path + ' \
                     -ref ' + reference_template_path + ' \
                     -out ' + template_flirtedRigid_path + ' \
                     -init ' + flirtRigid_path + ' \
                     -applyxfm \
                     -interp trilinear \
                     -verbose 1')

    # Apply rigid FLIRT to annotation
    print('FLIRT affine apply to annotation start')
    os.system('flirt -in ' + annotation_path + ' \
                     -ref ' + reference_template_path + ' \
                     -out ' + annotation_flirtedRigid_path + ' \
                     -init ' + flirtRigid_path + ' \
                     -applyxfm \
                     -interp nearestneighbour \
                     -verbose 1')

