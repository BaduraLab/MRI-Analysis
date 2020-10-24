from functions import imageFLIRT2defField
import nibabel as nib
from compress_pickle import dump, load
import numpy as np
import os
import glob
from scipy.stats import ttest_ind

# Define paths
data_path = os.path.join('Data', 'Mouse', 'Processed')
reference_path = os.path.join('Data', 'Mouse', 'Reference')
data_path_list = glob.glob(os.path.join(data_path, '*'))
ref_path = os.path.join(reference_path, 'annotation_50_reoriented_flirted_cropped.nii.gz')
ref_image = nib.load(ref_path)

subject_list = list()
genotype_list = list()
nData = len(data_path_list)
defField_magnitude_flirtRigid = np.empty(list(ref_image.shape)+[nData])
defField_magnitude_flirt = np.empty(list(ref_image.shape)+[nData])
defField_magnitude_syn = np.empty(list(ref_image.shape)+[nData])
for iData, Path in enumerate(data_path_list):
    print(iData)
    print(Path)

    subject = Path.split(os.sep)[-1]
    subject_list.append(subject)
    genotype_list.append(subject.split('_')[0])
    print(subject)

    # Define subject specific paths
    or_path = os.path.join(Path, subject + '.nii.gz')
    flirtRigid_path = os.path.join(Path, subject + '_flirtRigid.mat')
    flirt_path = os.path.join(Path, subject + '_flirt.mat')
    invflirt_path = os.path.join(Path, subject + '_invflirt.mat')
    flirtRigid_path = os.path.join(Path, subject + '_flirtRigid.mat')
    invflirtRigid_path = os.path.join(Path, subject + '_invflirtRigid.mat')
    syn_path = os.path.join(Path, subject + '_flirted_syn.pickle.gz')



    # rigid flirt vector field and magnitude
    defField_path = os.path.join(Path, subject + '_flirtRigid_defField.nii.gz')
    defField_magnitude_path = os.path.join(Path, subject + '_flirtRigid_defField_magnitude.nii.gz')
    with open(flirtRigid_path, 'r') as f:
        txt = f.read()
        flirtRigid = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
    print(flirtRigid)
    with open(invflirtRigid_path, 'r') as f:
        txt = f.read()
        invflirtRigid = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
    print(invflirtRigid)
    defField = imageFLIRT2defField(ref_image, invflirtRigid)
    defField_image = nib.Nifti1Image(defField, ref_image.affine)
    nib.save(defField_image, defField_path)
    defField_magnitude = np.sqrt(np.sum(np.power(defField, 2), axis=3))
    defField_magnitude_image = nib.Nifti1Image(defField_magnitude, ref_image.affine)
    nib.save(defField_magnitude_image, defField_magnitude_path)
    defField_magnitude_flirtRigid[:,:,:,iData] = defField_magnitude # Assign syn defField_magnitude's to 4D array


    # affine flirt vector field and magnitude
    defField_path = os.path.join(Path, subject + '_flirt_defField.nii.gz')
    defField_magnitude_path = os.path.join(Path, subject + '_flirt_defField_magnitude.nii.gz')
    with open(flirt_path, 'r') as f:
        txt = f.read()
        flirt = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
    print(flirt)
    with open(invflirt_path, 'r') as f:
        txt = f.read()
        invflirt = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
    print(invflirt)
    defField = imageFLIRT2defField(ref_image, invflirt)
    defField_image = nib.Nifti1Image(defField, ref_image.affine)
    nib.save(defField_image, defField_path)
    defField_magnitude = np.sqrt(np.sum(np.power(defField, 2), axis=3))
    defField_magnitude_image = nib.Nifti1Image(defField_magnitude, ref_image.affine)
    nib.save(defField_magnitude_image, defField_magnitude_path)
    defField_magnitude_flirt[:,:,:,iData] = defField_magnitude # Assign syn defField_magnitude's to 4D array



    # Calculate inverse flirt's and check whether it matches with read inverse flirt's



    # SyN vector field and magnitude
    defField_path = os.path.join(Path, subject + '_flirted_syn_defField.nii.gz')
    defField_magnitude_path = os.path.join(Path, subject + '_flirted_syn_defField_magnitude.nii.gz')

    # Load SyN of subject and get inverse field
    with open(syn_path, 'rb') as f:
        [mapping, metric, level_iters, sdr] = load(f, compression='gzip')
    defField = mapping.get_backward_field()
    del mapping, metric, level_iters, sdr

    defField_image = nib.Nifti1Image(defField, ref_image.affine)
    nib.save(defField_image, defField_path)
    defField_magnitude = np.sqrt(np.sum(np.power(defField, 2), axis=3))
    defField_magnitude_image = nib.Nifti1Image(defField_magnitude, ref_image.affine)
    nib.save(defField_magnitude_image, defField_magnitude_path)
    defField_magnitude_syn[:,:,:,iData] = defField_magnitude # Assign syn defFIeld_magnitude's to 4D array

# Average defField magnitude per genotype
genotype_WT_logical = np.array([genotype == 'WT' for genotype in genotype_list])
genotype_KO_logical = np.logical_not(genotype_WT_logical)
defField_magnitude_flirtRigid_meanWT = np.empty(defField_magnitude.shape[0:3])
defField_magnitude_flirt_meanWT = np.empty(defField_magnitude.shape[0:3])
defField_magnitude_syn_meanWT = np.empty(defField_magnitude.shape[0:3])
defField_magnitude_flirtRigid_meanKO = np.empty(defField_magnitude.shape[0:3])
defField_magnitude_flirt_meanKO = np.empty(defField_magnitude.shape[0:3])
defField_magnitude_syn_meanKO = np.empty(defField_magnitude.shape[0:3])
defField_magnitude_flirtRigid_pval = np.empty(defField_magnitude.shape[0:3])
defField_magnitude_flirt_pval = np.empty(defField_magnitude.shape[0:3])
defField_magnitude_syn_pval = np.empty(defField_magnitude.shape[0:3])
for i in range(defField_magnitude.shape[0]):
    for j in range(defField_magnitude.shape[1]):
        for k in range(defField_magnitude.shape[2]):
            defField_magnitude_flirtRigid_genotype = np.squeeze(defField_magnitude_flirtRigid[i,j,k,:])
            defField_magnitude_flirt_genotype = np.squeeze(defField_magnitude_flirt[i,j,k,:])
            defField_magnitude_syn_genotype = np.squeeze(defField_magnitude_syn[i,j,k,:])

            defField_magnitude_flirtRigid_meanWT[i,j,k] = np.mean(defField_magnitude_flirtRigid_genotype[genotype_WT_logical])
            defField_magnitude_flirt_meanWT[i,j,k] = np.mean(defField_magnitude_flirt_genotype[genotype_WT_logical])
            defField_magnitude_syn_meanWT[i,j,k]   = np.mean(defField_magnitude_syn_genotype[  genotype_WT_logical])
            defField_magnitude_flirtRigid_meanKO[i,j,k] = np.mean(defField_magnitude_flirtRigid_genotype[genotype_KO_logical])
            defField_magnitude_flirt_meanKO[i,j,k] = np.mean(defField_magnitude_flirt_genotype[genotype_KO_logical])
            defField_magnitude_syn_meanKO[i,j,k]   = np.mean(defField_magnitude_syn_genotype[  genotype_KO_logical])

            [_, defField_magnitude_flirtRigid_pval[i,j,k]] = ttest_ind(defField_magnitude_flirtRigid_genotype[genotype_WT_logical],
                                                                      defField_magnitude_flirtRigid_genotype[genotype_KO_logical],
                                                                      equal_var=False)
            [_, defField_magnitude_flirt_pval[i,j,k]] = ttest_ind(defField_magnitude_flirt_genotype[genotype_WT_logical],
                                                                  defField_magnitude_flirt_genotype[genotype_KO_logical],
                                                                  equal_var=False)
            [_, defField_magnitude_syn_pval[i,j,k]] = ttest_ind(defField_magnitude_syn_genotype[genotype_WT_logical],
                                                                defField_magnitude_syn_genotype[genotype_KO_logical],
                                                                equal_var=False)

# Save flirtRigid
defField_magnitude_meanWT_path = os.path.join(Path, subject + '_flirtRigid_defField_magnitude_meanWT.nii.gz')
defField_magnitude_meanKO_path = os.path.join(Path, subject + '_flirtRigid_defField_magnitude_meanKO.nii.gz')
defField_magnitude_pval_path = os.path.join(Path, subject + '_flirtRigid_defField_magnitude_pval.nii.gz')
defField_image = nib.Nifti1Image(defField_magnitude_flirtRigid_meanWT, ref_image.affine)
nib.save(defField_image, defField_magnitude_meanWT_path)
defField_image = nib.Nifti1Image(defField_magnitude_flirtRigid_meanKO, ref_image.affine)
nib.save(defField_image, defField_magnitude_meanKO_path)
defField_image = nib.Nifti1Image(defField_magnitude_flirtRigid_pval, ref_image.affine)
nib.save(defField_image, defField_magnitude_pval_path)

# Save flirt
defField_magnitude_meanWT_path = os.path.join(Path, subject + '_flirt_defField_magnitude_meanWT.nii.gz')
defField_magnitude_meanKO_path = os.path.join(Path, subject + '_flirt_defField_magnitude_meanKO.nii.gz')
defField_magnitude_pval_path = os.path.join(Path, subject + '_flirt_defField_magnitude_pval.nii.gz')
defField_image = nib.Nifti1Image(defField_magnitude_flirt_meanWT, ref_image.affine)
nib.save(defField_image, defField_magnitude_meanWT_path)
defField_image = nib.Nifti1Image(defField_magnitude_flirt_meanKO, ref_image.affine)
nib.save(defField_image, defField_magnitude_meanKO_path)
defField_image = nib.Nifti1Image(defField_magnitude_flirt_pval, ref_image.affine)
nib.save(defField_image, defField_magnitude_pval_path)

# Save syn
defField_magnitude_meanWT_path = os.path.join(Path, subject + '_syn_defField_magnitude_meanWT.nii.gz')
defField_magnitude_meanKO_path = os.path.join(Path, subject + '_syn_defField_magnitude_meanKO.nii.gz')
defField_magnitude_pval_path = os.path.join(Path, subject + '_syn_defField_magnitude_pval.nii.gz')
defField_image = nib.Nifti1Image(defField_magnitude_syn_meanWT, ref_image.affine)
nib.save(defField_image, defField_magnitude_meanWT_path)
defField_image = nib.Nifti1Image(defField_magnitude_syn_meanKO, ref_image.affine)
nib.save(defField_image, defField_magnitude_meanKO_path)
defField_image = nib.Nifti1Image(defField_magnitude_syn_pval, ref_image.affine)
nib.save(defField_image, defField_magnitude_pval_path)
