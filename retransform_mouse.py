# retransform mouse native annotation, flirtRigid, flirted
# unlike with human transform this would not have been necessary if synned annotation was saved
# for visualization purposes and consistency the retransform is sufficient and somewhat warranted respectively
# possibly need to ask additional mask
import os
import numpy as np
import glob
import datetime
from functions import zeroPadImage
import nibabel as nib
from functions import save_image
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.imwarp import DiffeomorphicMap
from dipy.align.metrics import CCMetric
from compress_pickle import dump, load

# Define
data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
subject_path_list = glob.glob(os.path.join(data_path, '*'))
reference_path = os.path.join('Data', 'Mouse', 'Reference')
reference_template_path = os.path.join(reference_path, 'average_template_25_reoriented_flirted.nii.gz')
reference_template_image = nib.load(reference_template_path)
reference_template = reference_template_image.get_fdata()

subject_path_list = [subject_path_list[0]]
for iSubjectPath, SubjectPath in enumerate(subject_path_list):
    print(iSubjectPath)
    print(iSubjectPath)
    print(SubjectPath)
    print(datetime.datetime.now())

    # Define subject template paths
    subject_string = SubjectPath.split(os.sep)[-1]
    template_path = os.path.join(SubjectPath, 'FLIRT', subject_string + '_inmasked.nii.gz')

    # Define subject native annotation path
    annotation_path_list = [os.path.join(SubjectPath, 'allen_annotation_invsynned_to_' + subject_string + '.nii.gz')]

    # create "retransform" folder
    retransform_folder = os.path.join(SubjectPath, 'retransform')
    if not os.path.exists(retransform_folder):
        os.makedirs(retransform_folder)
    print(retransform_folder)

    template_flirtedRigid_path = os.path.join(retransform_folder, subject_string + '_flirtedRigid.nii.gz')
    flirtRigid_path = os.path.join(retransform_folder, 'flirtRigid.mat')
    template_flirted_path = os.path.join(retransform_folder, subject_string + '_flirted.nii.gz')
    flirt_path = os.path.join(retransform_folder, 'flirt.mat')
    template_synned_path = os.path.join(retransform_folder, subject_string + '_synned.nii.gz')
    syn_path = os.path.join(retransform_folder, 'syn.pickle.gz')



    ## Compute rigid flirt, flirt and SyN by transforming templates to the refernce template
    # Rigid FLIRT inmasked subject to reference
    print('FLIRT rigid start')
    os.system('flirt -in ' + template_path + ' \
                     -ref ' + reference_template_path + ' \
                     -out ' + template_flirtedRigid_path + ' \
                     -omat ' + flirtRigid_path + ' \
                     -dof ' + '6' + ' \
                     -verbose 0')    # FLIRT subject to reference
    #                  # -dof ' + '6' + ' \

    # Affine FLIRT inmasked subject to reference
    print('FLIRT affine start')
    os.system('flirt -in ' + template_flirtedRigid_path + ' \
                     -ref ' + reference_template_path + ' \
                     -out ' + template_flirted_path + ' \
                     -omat ' + flirt_path + ' \
                     -verbose 0')    # FLIRT subject to reference
    template_flirted_image = nib.load(template_flirted_path)
    template_flirted = template_flirted_image.get_fdata()

    ## SyN
    print('SyN start')
    metric = CCMetric(3)
    level_iters = [10, 10, 5, 5, 5]
    sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)
    mapping = sdr.optimize(static=reference_template,
                           moving=template_flirted,
                           static_grid2world=reference_template_image.get_qform(),
                           moving_grid2world=template_flirted_image.get_qform())
    with open(syn_path, 'wb') as f:
        dump([mapping, metric, level_iters, sdr], f, protocol=4, compression='gzip')
    template_synned = mapping.transform(template_flirted)
    template_synned_image = nib.Nifti1Image(template_synned, np.eye(4))
    template_synned_image.set_qform(template_flirted_image.get_qform(), code=1)
    template_synned_image.set_sform(np.eye(4), code=0)
    nib.save(template_synned_image, template_synned_path)



    ## Apply rigid flirt, flirt and SyN to annotations
    print('Apply transformations to annotations start')
    for annotation_path in annotation_path_list:

        # flirt Rigid
        annotation_flirtedRigid_path = os.path.join(retransform_folder, annotation_path.split(os.sep)[-1].split('.')[0] + '_flirtedRigid.nii.gz')
        os.system('flirt -in ' + annotation_path + ' \
                         -ref ' + reference_template_path + ' \
                         -out ' + annotation_flirtedRigid_path + ' \
                         -init ' + flirtRigid_path + ' \
                         -applyxfm \
                         -interp ' + 'nearestneighbour' + ' \
                         -verbose 0')

        # flirt
        annotation_flirted_path = os.path.join(retransform_folder, annotation_path.split(os.sep)[-1].split('.')[0] + '_flirted.nii.gz')
        os.system('flirt -in ' + annotation_flirtedRigid_path + ' \
                         -ref ' + reference_template_path + ' \
                         -out ' + annotation_flirted_path + ' \
                         -init ' + flirt_path + ' \
                         -applyxfm \
                         -interp ' + 'nearestneighbour' + ' \
                         -verbose 0')
        annotation_flirted_image = nib.load(annotation_flirted_path)
        annotation_flirted = annotation_flirted_image.get_fdata()

        # SyN
        annotation_synned_path = os.path.join(retransform_folder, annotation_path.split(os.sep)[-1].split('.')[0] + '_synned.nii.gz')
        annotation_synned = mapping.transform(annotation_flirted, interpolation='nearest')
        annotation_synned_image = nib.Nifti1Image(annotation_synned, np.eye(4))
        annotation_synned_image.set_qform(annotation_flirted_image.get_qform(), code=1)
        annotation_synned_image.set_sform(np.eye(4), code=0)
        nib.save(annotation_synned_image, annotation_synned_path)
