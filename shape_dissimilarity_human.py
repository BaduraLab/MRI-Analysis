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

# Define paths
data_path = os.path.join('Data', 'Human', 'Processed')
reference_path = os.path.join('Data', 'Human', 'Reference')
analysis_path = os.path.join('Data', 'Human', 'Analysis')
data_path_list = glob.glob(os.path.join(data_path, '*'))
reference_annotation_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT.nii.gz'),
                        os.path.join(reference_path, 'subcortical', 'prob_atlas_bilateral_thrarg_0.4.nii.gz')]
reference_template_path_list = [os.path.join(reference_path, 'suit', 'templates', 'SUIT.nii'),
                      os.path.join(reference_path, 'subcortical', 'CIT168_T1w_700um_reoriented.nii.gz')]
reference_structure_path_list = [os.path.join(reference_path, 'suit', 'atlasesSUIT', 'Lobules-SUIT_mc.csv'),
                       os.path.join(reference_path, 'subcortical', 'subcortical.csv')]
reference_structure_list = [pd.read_csv(reference_structure_path_list[0]),
                            pd.read_csv(reference_structure_path_list[1])]
id_list_list = [[1], [7]]
reference_annotation_name_list = ['orsuit',
                                  'subcortical']

RefIsoPath = os.path.join(reference_path, 'Isolations')
if not os.path.exists(RefIsoPath):
    os.makedirs(RefIsoPath)
zeroPadRatio = 0.5

# loop over atlases
# use name instead of acronym
# get volumeIntegers differently compared to id_custom
# all paths need to be checked



for iGroup in range(len(reference_annotation_path_list)):
    print(f'iGroup = {iGroup}')
    reference_annotation_path = reference_annotation_path_list[iGroup]
    reference_template_path = reference_template_path_list[iGroup]
    reference_structure = reference_structure_list[iGroup]

    reference_template_image = nib.load(reference_template_path)
    reference_template = reference_template_image.get_fdata()
    # reference_shape = reference_template.shape
    reference_annotation_image = nib.load(reference_annotation_path)
    reference_annotation = reference_annotation_image.get_fdata()
    reference_annotation = np.round(reference_annotation).astype(int)

    id_list = id_list_list[iGroup]
    for id in id_list:
        structure_name = reference_structure.loc[reference_structure['id_custom'] == id, 'name'].iloc[0]
        start_time = datetime.datetime.now()
        print(f'Start time = {start_time}')

        # isolate reference template and annotation
        reference_template_zeropadded_path = os.path.join(RefIsoPath, 'reference_template_' + reference_annotation_name_list[iGroup] + '_' + structure_name + '_zeropadded.nii.gz')
        reference_annotation_in_group = np.isin(reference_annotation, id)
        reference_template_group = reference_template * reference_annotation_in_group
        reference_template_zeropadded, reference_crop_index = zeroPadImage(reference_template_group, reference_annotation_in_group, zeroPadRatio)
        save_image(reference_template_zeropadded, reference_template_image, reference_template_zeropadded_path)
        reference_shape = reference_template_zeropadded.shape

        # Define subject related paths
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
            template_zeropadded_path = os.path.join(IsoPath, subject + '_inmasked_' + structure_name + '_zeropadded.nii.gz')
            flirtRigid_path = os.path.join(IsoPath, subject + '_' + structure_name + '_flirtRigidSD.mat')
            flirtedRigid_path = os.path.join(IsoPath, subject + '_' + structure_name + '_flirtedRigidSD.nii.gz')
            flirtAffine_path = os.path.join(IsoPath, subject + '_' + structure_name + '_flirtAffineSD.mat')
            flirtedAffine_path = os.path.join(IsoPath, subject + '_' + structure_name + '_flirtedAffineSD.nii.gz')
            invflirtRigid_path = os.path.join(IsoPath, subject + '_' + structure_name + '_invflirtRigidSD.mat')
            invflirtAffine_path = os.path.join(IsoPath, subject + '_' + structure_name + '_invflirtAffineSD.mat')
            syn_path = os.path.join(IsoPath, subject + '_' + structure_name + '_synSD.pickle.gz') # also include crop_index, flirtRigid and flirtAffine
            synned_path = os.path.join(IsoPath, subject + '_' + structure_name + '_synnedSD.nii.gz')
            jacobian_det_flirtRigid_path = os.path.join(IsoPath, subject + '_' + structure_name + '_flirtRigidSD.nii.gz')
            jacobian_det_flirtAffine_path = os.path.join(IsoPath, subject + '_' + structure_name + '_flirtAffineSD.nii.gz')
            jacobian_det_syn_path = os.path.join(IsoPath, subject + '_' + structure_name + '_synSD.nii.gz')
            jacobian_det_all_path = os.path.join(IsoPath, subject + '_' + structure_name + '_allSD.nii.gz')
            jacobian_det_VC_path = os.path.join(IsoPath, subject + '_' + structure_name + '_VCSD.nii.gz')
            jacobian_det_flirt_path = os.path.join(IsoPath, subject + '_' + structure_name + '_flirtSD.nii.gz')

            # Load native template and annotation
            template_path = glob.glob(os.path.join(DataPath, subject + '_reoriented.nii.gz'))
            template_path = template_path[0]
            template_image = nib.load(template_path)
            template = template_image.get_fdata()
            print(f'template_path = {template_path}')
            if reference_annotation_name_list[iGroup] == 'orsuit':
                annotation_path = os.path.join(DataPath, subject + '_annotation_orsuit_thrarg_adjusted_lobular.nii.gz')
            elif reference_annotation_name_list[iGroup] == 'subcortical':
                annotation_path = os.path.join(DataPath, subject + '_annotation_subcortical_thrarg.nii.gz')
            annotation_image = nib.load(annotation_path)
            annotation = annotation_image.get_fdata()
            annotation = np.round(annotation).astype(int)
            print(f'annotation_path = {annotation_path}')

            # isolate structures, zeropad native template, reference template
            annotation_in_group = np.isin(annotation, id)
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
            defField_flirtRigid = imageFLIRT2defField(reference_template_image, invflirtRigid)

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
            defField_flirtAffine = imageFLIRT2defField(reference_template_image, invflirtAffine)

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

            # Calculate jacobian determinants
            jacobian_det_flirtRigid = computeJacobianDet(defField_flirtRigid)
            jacobian_det_flirtAffine = computeJacobianDet(defField_flirtAffine)
            jacobian_det_syn = computeJacobianDet(back_field)
            jacobian_det_all = computeJacobianDet(defField_flirtRigid + defField_flirtAffine + back_field)
            jacobian_det_VC = computeJacobianDet(defField_flirtAffine + back_field)
            jacobian_det_flirt = computeJacobianDet(defField_flirtRigid + defField_flirtAffine)

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

            print(f'Processing time = {datetime.datetime.now() - start_time}')



        # Average defField magnitude per genotype
        start_time = datetime.datetime.now()
        print(f'Start time for post genotype averaging = {start_time}')
        genotype_WT_logical = np.array([genotype == 'WT' for genotype in genotype_list])
        genotype_KO_logical = np.logical_not(genotype_WT_logical)
        defField_Jdet_flirtRigid_meanWT = np.empty(reference_shape)
        defField_Jdet_flirtAffine_meanWT = np.empty(reference_shape)
        defField_Jdet_syn_meanWT = np.empty(reference_shape)
        defField_Jdet_flirtRigid_meanKO = np.empty(reference_shape)
        defField_Jdet_flirtAffine_meanKO = np.empty(reference_shape)
        defField_Jdet_syn_meanKO = np.empty(reference_shape)
        defField_Jdet_flirtRigid_pval = np.empty(reference_shape)
        defField_Jdet_flirtAffine_pval = np.empty(reference_shape)
        defField_Jdet_syn_pval = np.empty(reference_shape)
        for i in range(reference_shape[0]):
            for j in range(reference_shape[1]):
                for k in range(reference_shape[2]):
                    defField_Jdet_flirtRigid_genotype = np.squeeze(defField_Jdet_flirtRigid[i, j, k, :])
                    defField_Jdet_flirtAffine_genotype = np.squeeze(defField_Jdet_flirtAffine[i, j, k, :])
                    defField_Jdet_syn_genotype = np.squeeze(defField_Jdet_syn[i, j, k, :])

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

            print(f'Processing time for post genotype averaging = {datetime.datetime.now() - start_time}')
