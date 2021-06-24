from functions import imageFLIRT2defField
import nibabel as nib
from compress_pickle import dump, load
import numpy as np
import os
import glob
from scipy.stats import ttest_ind
import numdifftools
from compress_pickle import dump, load
from dipy.align.metrics import CCMetric
import datetime
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
import sympy
from functions import computeJacobianDet
from functions import save_image
from functions import unzeroPadImage
from functions import zeroPadImage
from functions import getChildStructures_mouse
import pandas as pd

## MAIN: we include manual annotation by isolation and superposition of Jacobian determinant arrays
# flirtRigid flirt SyN manual (native-reference)
# isolate manually adjusted structure in reference space
# flirtRigid flirt SyN (reference-native)
# 3 vector fields
# can be combined in different ways (7, most interesting one is all put together and affine and non-linear put together --> volume changing)
# calculate jacobian vector field --> matrix field
# calculate determinant
# only save determinant scalar field

## SECOND: we do not include manual annotation and simply use the inverse transforms that were actually applied as vector fields
# problem is that no flirtRigid used originally, only affine and non-linear

## For each mouse and each structure you would need to do flirtRigid, flirtAffine and SyN as well as 3 Jacobian determinant calculations which loops through each voxel
# Then possibly you can do an additional 3 Jacobian determinant calculations with "all", "VC" and "flirt" combinations.
# In the end the jacobian determinants can be saved into a 4D array with 4th dimension being mice, and then averaged according to genotype.
# Cropping of areas will increase computational efficiency of Jacobian step, but possibly not of flirt and SyN steps

# Define paths
data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
reference_path = os.path.join('Data', 'Mouse', 'Reference')
analysis_path = os.path.join('Data', 'Mouse', 'Analysis')
data_path_list = glob.glob(os.path.join(data_path, '*'))
# data_path_list = [data_path_list[0]]
reference_annotation_path = os.path.join(reference_path, 'annotation_25_reoriented_flirted.nii.gz')
reference_annotation_image = nib.load(reference_annotation_path)
reference_annotation = reference_annotation_image.get_fdata()
reference_annotation = reference_annotation = np.round(reference_annotation).astype(int)
reference_template_path = os.path.join(reference_path, 'average_template_25_reoriented_flirted.nii.gz')
reference_template_image = nib.load(reference_template_path)
reference_template = reference_template_image.get_fdata()
# reference_shape = reference_template.shape
reference_structure_path = os.path.join(reference_path, 'structure_graph_plus.csv')
reference_structure = pd.read_csv(reference_structure_path)

RefIsoPath = os.path.join(reference_path, 'Isolations')
if not os.path.exists(RefIsoPath):
    os.makedirs(RefIsoPath)
zeroPadRatio = 0.3 # 0.3 worked for SyN, make it 0.5 with some margin of error

id_list = list(np.unique(reference_structure['id_custom']))
id_list = [2001]
id_list_list = [[54, 268, 418]]
acronym_list = ['SN']
for iId, id_list in enumerate(id_list_list):
# for id in id_list:
    # structure_acronym = reference_structure.loc[reference_structure['id_custom'] == id, 'acronym'].iloc[0]
    structure_acronym = acronym_list[iId]
    start_time_id = datetime.datetime.now()
    print(f'Start time for {structure_acronym} = {start_time_id}')

    # isolate structure, zeropad reference template and annotation
    reference_template_zeropadded_path = os.path.join(RefIsoPath, 'average_template_25_reoriented_flirted_' + structure_acronym + '_zeropadded.nii.gz')
    # child = getChildStructures_mouse(id, reference_structure)
    # reference_annotation_in_group = np.isin(reference_annotation, child)
    reference_annotation_in_group = np.isin(reference_annotation, id_list)
    reference_template_group = reference_template * reference_annotation_in_group
    reference_template_zeropadded, reference_crop_index = zeroPadImage(reference_template_group,
                                                                       reference_annotation_in_group, zeroPadRatio)
    reference_template_zeropadded_image = save_image(reference_template_zeropadded, reference_template_image,
                                                     reference_template_zeropadded_path)
    reference_shape = reference_template_zeropadded.shape

    # Define subject related lists
    subject_list = list()
    genotype_list = list()
    nData = len(data_path_list)
    defField_Jdet_flirtRigid = np.empty(list(reference_shape) + [nData])
    defField_Jdet_flirtAffine = np.empty(list(reference_shape) + [nData])
    defField_Jdet_syn = np.empty(list(reference_shape) + [nData])
    defField_Jdet_all = np.empty(list(reference_shape) + [nData])
    defField_Jdet_VC = np.empty(list(reference_shape) + [nData])
    defField_Jdet_flirt = np.empty(list(reference_shape) + [nData])
    for iDataPath, DataPath in enumerate(data_path_list):
        print(f'iDataPath = {iDataPath}')
        print(f'DataPath = {DataPath}')

        # Define subject specific paths
        IsoPath = os.path.join(DataPath, 'Isolations', 'SD')
        if not os.path.exists(IsoPath):
            os.makedirs(IsoPath)
        subject = DataPath.split(os.sep)[-1]
        subject_list.append(subject)
        genotype = subject.split('_')[0]
        genotype_list.append(genotype)
        print(f'subject = {subject}')
        print(f'genotype = {genotype}')
        start_time = datetime.datetime.now()
        print(f'Start time for {subject} processing {structure_acronym} = {start_time}')
        template_zeropadded_path = os.path.join(IsoPath, subject + '_inmasked_' + structure_acronym + '_zeropadded.nii.gz')
        flirtRigid_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_flirtRigidSD.mat')
        flirtedRigid_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_flirtedRigidSD.nii.gz')
        flirtAffine_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_flirtAffineSD.mat')
        flirtedAffine_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_flirtedAffineSD.nii.gz')
        invflirtRigid_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_invflirtRigidSD.mat')
        invflirtAffine_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_invflirtAffineSD.mat')
        syn_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_synSD.pickle.gz') # also include crop_index, flirtRigid and flirtAffine
        synned_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_synnedSD.nii.gz')
        jacobian_det_flirtRigid_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_flirtRigidSD.nii.gz')
        jacobian_det_flirtAffine_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_flirtAffineSD.nii.gz')
        jacobian_det_syn_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_synSD.nii.gz')
        jacobian_det_all_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_allSD.nii.gz')
        jacobian_det_VC_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_VCSD.nii.gz')
        jacobian_det_flirt_path = os.path.join(IsoPath, subject + '_' + structure_acronym + '_flirtSD.nii.gz')

        # Load native template and annotation
        template_path = os.path.join(DataPath, 'FLIRT', subject + '_inmasked' + '.nii.gz')
        template_image = nib.load(template_path)
        template = template_image.get_fdata()
        print(f'template_path = {template_path}')
        annotation_path = os.path.join(DataPath, 'allen_annotation_invsynned_to_' + subject + '_adjusted_cerebellum_lobular.nii.gz')
        annotation_image = nib.load(annotation_path)
        annotation = annotation_image.get_fdata()
        annotation = np.round(annotation).astype(int)
        print(f'annotation_path = {annotation_path}')

        # isolate structures, zeropad native template, reference template
        annotation_in_group = np.isin(annotation, id_list)
        template_group = template * annotation_in_group
        template_zeropadded, native_crop_index = zeroPadImage(template_group, annotation_in_group, zeroPadRatio)
        save_image(template_zeropadded, template_image, template_zeropadded_path)

        # flirtRigid template_path to reference_path
        print('FLIRT rigid start')
        os.system('flirt -in ' + template_zeropadded_path + ' \
                         -ref ' + reference_template_zeropadded_path + ' \
                         -out ' + flirtedRigid_path + ' \
                         -omat ' + flirtRigid_path + ' \
                         -dof ' + '6' + ' \
                         -verbose 0')  # FLIRT subject to reference
        os.system('convert_xfm -omat ' + invflirtRigid_path + ' -inverse ' + flirtRigid_path)
        with open(flirtRigid_path, 'r') as f:
            txt = f.read()
            flirtRigid = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
        print(flirtRigid)
        with open(invflirtRigid_path, 'r') as f:
            txt = f.read()
            invflirtRigid = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
        print(f'invflirtRigid = {invflirtRigid}')
        defField_flirtRigid = imageFLIRT2defField(reference_template_zeropadded_image, invflirtRigid)

        # flirtAffine template_flirtedRigid_path to reference_path
        print('FLIRT affine start')
        os.system('flirt -in ' + flirtedRigid_path + ' \
                         -ref ' + reference_template_zeropadded_path + ' \
                         -out ' + flirtedAffine_path + ' \
                         -omat ' + flirtAffine_path + ' \
                         -verbose 0')
        flirtedAffine_image = nib.load(flirtedAffine_path)
        flirtedAffine = flirtedAffine_image.get_fdata()
        os.system('convert_xfm -omat ' + invflirtAffine_path + ' -inverse ' + flirtAffine_path)
        with open(flirtAffine_path, 'r') as f:
            txt = f.read()
            flirtAffine = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
        print(flirtAffine)
        with open(invflirtAffine_path, 'r') as f:
            txt = f.read()
            invflirtAffine = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
        print(f'invflirtAffine = {invflirtAffine}')
        defField_flirtAffine = imageFLIRT2defField(reference_template_zeropadded_image, invflirtAffine)

        # SyN template_flirtedAffine_path to reference_path
        print('SyN')
        metric = CCMetric(3)
        level_iters = [10, 10, 5, 5, 5]
        sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)
        mapping = sdr.optimize(static=reference_template_zeropadded,
                               moving=flirtedAffine,
                               static_grid2world=reference_template_image.get_qform(),
                               moving_grid2world=flirtedAffine_image.get_qform())
        with open(syn_path, 'wb') as f:
            dump([mapping, metric, level_iters, sdr,
                  native_crop_index, reference_crop_index,
                  flirtRigid, invflirtRigid, flirtAffine, invflirtAffine],
                 f, protocol=4, compression='gzip')
        # with open(mouse_masked_syn_path, 'rb') as f:
        #     [mapping, metric, level_iters, sdr] = load(f)
        # Check SyN
        forw_field = mapping.get_forward_field()
        back_field = mapping.get_backward_field()
        forw_SS = np.sum(np.power(forw_field, 2))
        back_SS = np.sum(np.power(back_field, 2))
        dif_SSD = np.sum(np.power(forw_field + back_field, 2))
        dif_SSD_norm = dif_SSD / ((forw_SS + back_SS) / 2)
        print(f'dif_SSD_norm = {dif_SSD_norm}')

        # JD calculation, saving images per subject and assignment of JD deformation fields
        start_time_JD = datetime.datetime.now()
        # Calculate jacobian determinants
        jacobian_det_flirtRigid = computeJacobianDet(defField_flirtRigid, reference_annotation_in_group)
        jacobian_det_flirtAffine = computeJacobianDet(defField_flirtAffine, reference_annotation_in_group)
        jacobian_det_syn = computeJacobianDet(back_field, reference_annotation_in_group)
        jacobian_det_all = computeJacobianDet(defField_flirtRigid + defField_flirtAffine + back_field, reference_annotation_in_group)
        jacobian_det_VC = computeJacobianDet(defField_flirtAffine + back_field, reference_annotation_in_group)
        jacobian_det_flirt = computeJacobianDet(defField_flirtRigid + defField_flirtAffine, reference_annotation_in_group)

        # Make images for all jacobian determinants
        save_image(jacobian_det_flirtRigid, reference_template_image, jacobian_det_flirtRigid_path)
        save_image(jacobian_det_flirtAffine, reference_template_image, jacobian_det_flirtAffine_path)
        save_image(jacobian_det_syn, reference_template_image, jacobian_det_syn_path)
        save_image(jacobian_det_all, reference_template_image, jacobian_det_all_path)
        save_image(jacobian_det_VC, reference_template_image, jacobian_det_VC_path)
        save_image(jacobian_det_flirt, reference_template_image, jacobian_det_flirt_path)

        # Assign Jacobian determinants
        defField_Jdet_flirtRigid[:, :, :, iDataPath] = jacobian_det_flirtRigid
        defField_Jdet_flirtAffine[:, :, :, iDataPath] = jacobian_det_flirtAffine
        defField_Jdet_syn[:, :, :, iDataPath] = jacobian_det_syn
        defField_Jdet_all[:, :, :, iDataPath] = jacobian_det_all
        defField_Jdet_VC[:, :, :, iDataPath] = jacobian_det_VC
        defField_Jdet_flirt[:, :, :, iDataPath] = jacobian_det_flirt

        print(f'JD processing time for {subject} processing {structure_acronym} = {datetime.datetime.now() - start_time_JD}')
        print(f'Processing time for {subject} processing {structure_acronym} = {datetime.datetime.now() - start_time}')



    # Average defField magnitude per genotype
    start_time = datetime.datetime.now()
    print(f'Start time for post genotype averaging = {start_time}')
    genotype_WT_logical = np.array([genotype == 'WT' for genotype in genotype_list])
    genotype_KO_logical = np.logical_not(genotype_WT_logical)
    defField_Jdet_flirtRigid_meanWT = np.zeros(reference_shape)
    defField_Jdet_flirtAffine_meanWT = np.zeros(reference_shape)
    defField_Jdet_syn_meanWT = np.zeros(reference_shape)
    defField_Jdet_flirtRigid_meanKO = np.zeros(reference_shape)
    defField_Jdet_flirtAffine_meanKO = np.zeros(reference_shape)
    defField_Jdet_syn_meanKO = np.zeros(reference_shape)
    defField_Jdet_flirtRigid_pval = np.zeros(reference_shape)
    defField_Jdet_flirtAffine_pval = np.zeros(reference_shape)
    defField_Jdet_syn_pval = np.zeros(reference_shape)
    for i in range(reference_shape[0]):
        for j in range(reference_shape[1]):
            for k in range(reference_shape[2]):
                if reference_annotation_in_group[i,j,k]:
                    defField_Jdet_flirtRigid_genotype = np.squeeze(defField_Jdet_flirtRigid[i, j, k, :])
                    defField_Jdet_flirtAffine_genotype = np.squeeze(defField_Jdet_flirtAffine[i, j, k, :])
                    defField_Jdet_syn_genotype = np.squeeze(defField_Jdet_syn[i, j, k, :])
                    defField_Jdet_flirtRigid_genotype = np.squeeze(defField_Jdet_flirtRigid[i, j, k, :])
                    defField_Jdet_flirtAffine_genotype = np.squeeze(defField_Jdet_flirtAffine[i, j, k, :])
                    defField_Jdet_syn_genotype = np.squeeze(defField_Jdet_syn[i, j, k, :])

                    defField_Jdet_flirtRigid_meanWT[i,j,k] = np.mean(defField_Jdet_flirtRigid_genotype[genotype_WT_logical])
                    defField_Jdet_flirtAffine_meanWT[i,j,k] = np.mean(defField_Jdet_flirtAffine_genotype[genotype_WT_logical])
                    defField_Jdet_syn_meanWT[i,j,k]   = np.mean(defField_Jdet_syn_genotype[  genotype_WT_logical])
                    defField_Jdet_flirtRigid_meanKO[i,j,k] = np.mean(defField_Jdet_flirtRigid_genotype[genotype_KO_logical])
                    defField_Jdet_flirtAffine_meanKO[i,j,k] = np.mean(defField_Jdet_flirtAffine_genotype[genotype_KO_logical])
                    defField_Jdet_syn_meanKO[i,j,k]   = np.mean(defField_Jdet_syn_genotype[  genotype_KO_logical])
                    defField_Jdet_flirtRigid_meanWT[i,j,k] = np.mean(defField_Jdet_flirtRigid_genotype[genotype_WT_logical])
                    defField_Jdet_flirtAffine_meanWT[i,j,k] = np.mean(defField_Jdet_flirtAffine_genotype[genotype_WT_logical])
                    defField_Jdet_syn_meanWT[i,j,k]   = np.mean(defField_Jdet_syn_genotype[  genotype_WT_logical])
                    defField_Jdet_flirtRigid_meanKO[i,j,k] = np.mean(defField_Jdet_flirtRigid_genotype[genotype_KO_logical])
                    defField_Jdet_flirtAffine_meanKO[i,j,k] = np.mean(defField_Jdet_flirtAffine_genotype[genotype_KO_logical])
                    defField_Jdet_syn_meanKO[i,j,k]   = np.mean(defField_Jdet_syn_genotype[  genotype_KO_logical])

                    [_, defField_Jdet_flirtRigid_pval[i,j,k]] = ttest_ind(defField_Jdet_flirtRigid_genotype[genotype_WT_logical],
                                                                               defField_Jdet_flirtRigid_genotype[genotype_KO_logical],
                                                                               equal_var=False)
                    [_, defField_Jdet_flirtAffine_pval[i,j,k]] = ttest_ind(defField_Jdet_flirtAffine_genotype[genotype_WT_logical],
                                                                          defField_Jdet_flirtAffine_genotype[genotype_KO_logical],
                                                                          equal_var=False)
                    [_, defField_Jdet_syn_pval[i,j,k]] = ttest_ind(defField_Jdet_syn_genotype[genotype_WT_logical],
                                                                        defField_Jdet_syn_genotype[genotype_KO_logical],
                                                                        equal_var=False)
                    [_, defField_Jdet_flirtRigid_pval[i,j,k]] = ttest_ind(defField_Jdet_flirtRigid_genotype[genotype_WT_logical],
                                                                               defField_Jdet_flirtRigid_genotype[genotype_KO_logical],
                                                                               equal_var=False)
                    [_, defField_Jdet_flirtAffine_pval[i,j,k]] = ttest_ind(defField_Jdet_flirtAffine_genotype[genotype_WT_logical],
                                                                          defField_Jdet_flirtAffine_genotype[genotype_KO_logical],
                                                                          equal_var=False)
                    [_, defField_Jdet_syn_pval[i,j,k]] = ttest_ind(defField_Jdet_syn_genotype[genotype_WT_logical],
                                                                        defField_Jdet_syn_genotype[genotype_KO_logical],
                                                                        equal_var=False)

    # Save flirtRigid
    defField_magnitude_meanWT_path = os.path.join(analysis_path, 'SD_flirtRigid_defField_Jdet_meanWT.nii.gz')
    defField_magnitude_meanKO_path = os.path.join(analysis_path, 'SD_flirtRigid_defField_Jdet_meanKO.nii.gz')
    defField_magnitude_pval_path = os.path.join(analysis_path, 'SD_flirtRigid_defField_Jdet_pval.nii.gz')
    defField_image = nib.Nifti1Image(defField_Jdet_flirtRigid_meanWT, reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_meanWT_path)
    defField_image = nib.Nifti1Image(defField_Jdet_flirtRigid_meanKO, reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_meanKO_path)
    defField_image = nib.Nifti1Image(np.abs(np.log10(defField_Jdet_flirtRigid_pval)), reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_pval_path)

    # Save flirt
    defField_magnitude_meanWT_path = os.path.join(analysis_path, 'SD_flirtAffine_defField_Jdet_meanWT.nii.gz')
    defField_magnitude_meanKO_path = os.path.join(analysis_path, 'SD_flirtAffine_defField_Jdet_meanKO.nii.gz')
    defField_magnitude_pval_path = os.path.join(analysis_path, 'SD_flirtAffine_defField_Jdet_pval.nii.gz')
    defField_image = nib.Nifti1Image(defField_Jdet_flirtAffine_meanWT, reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_meanWT_path)
    defField_image = nib.Nifti1Image(defField_Jdet_flirtAffine_meanKO, reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_meanKO_path)
    defField_image = nib.Nifti1Image(np.abs(np.log10(defField_Jdet_flirtAffine_pval)), reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_pval_path)

    # Save syn
    defField_magnitude_meanWT_path = os.path.join(analysis_path, 'SD_syn_defField_magnitude_meanWT.nii.gz')
    defField_magnitude_meanKO_path = os.path.join(analysis_path, 'SD_syn_defField_magnitude_meanKO.nii.gz')
    defField_magnitude_pval_path = os.path.join(analysis_path, 'SD_syn_defField_magnitude_pval.nii.gz')
    defField_image = nib.Nifti1Image(defField_Jdet_syn_meanWT, reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_meanWT_path)
    defField_image = nib.Nifti1Image(defField_Jdet_syn_meanKO, reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_meanKO_path)
    defField_image = nib.Nifti1Image(np.abs(np.log10(defField_Jdet_syn_pval)), reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_pval_path)

    # Save syn
    defField_magnitude_meanWT_path = os.path.join(analysis_path, 'SD_syn_defField_magnitude_meanWT.nii.gz')
    defField_magnitude_meanKO_path = os.path.join(analysis_path, 'SD_syn_defField_magnitude_meanKO.nii.gz')
    defField_magnitude_pval_path = os.path.join(analysis_path, 'SD_syn_defField_magnitude_pval.nii.gz')
    defField_image = nib.Nifti1Image(defField_Jdet_syn_meanWT, reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_meanWT_path)
    defField_image = nib.Nifti1Image(defField_Jdet_syn_meanKO, reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_meanKO_path)
    defField_image = nib.Nifti1Image(np.abs(np.log10(defField_Jdet_syn_pval)), reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_pval_path)

    # Save syn
    defField_magnitude_meanWT_path = os.path.join(analysis_path, 'SD_syn_defField_magnitude_meanWT.nii.gz')
    defField_magnitude_meanKO_path = os.path.join(analysis_path, 'SD_syn_defField_magnitude_meanKO.nii.gz')
    defField_magnitude_pval_path = os.path.join(analysis_path, 'SD_syn_defField_magnitude_pval.nii.gz')
    defField_image = nib.Nifti1Image(defField_Jdet_syn_meanWT, reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_meanWT_path)
    defField_image = nib.Nifti1Image(defField_Jdet_syn_meanKO, reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_meanKO_path)
    defField_image = nib.Nifti1Image(np.abs(np.log10(defField_Jdet_syn_pval)), reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_pval_path)

    # Save syn
    defField_magnitude_meanWT_path = os.path.join(analysis_path, 'SD_syn_defField_magnitude_meanWT.nii.gz')
    defField_magnitude_meanKO_path = os.path.join(analysis_path, 'SD_syn_defField_magnitude_meanKO.nii.gz')
    defField_magnitude_pval_path = os.path.join(analysis_path, 'SD_syn_defField_magnitude_pval.nii.gz')
    defField_image = nib.Nifti1Image(defField_Jdet_syn_meanWT, reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_meanWT_path)
    defField_image = nib.Nifti1Image(defField_Jdet_syn_meanKO, reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_meanKO_path)
    defField_image = nib.Nifti1Image(np.abs(np.log10(defField_Jdet_syn_pval)), reference_template_image.affine)
    nib.save(defField_image, defField_magnitude_pval_path)

    print(f'Processing time for post genotype averaging = {datetime.datetime.now() - start_time}')
    print(f'Processing time for {structure_acronym} = {datetime.datetime.now() - start_time_id}')
