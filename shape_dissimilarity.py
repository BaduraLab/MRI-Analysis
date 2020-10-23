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
ref_path = os.path.join(reference_path, 'annotation_25_reoriented_flirted.nii.gz')
ref_image = nib.load(ref_path)

subject_list = list()
genotype_list = list()
nData = len(data_path_list)
defField_magnitude_flirt = np.empty(shape=ref_image.shape+[nData])
defField_magnitude_syn = np.empty(shape=ref_image.shape+[nData])
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
    defField_path = os.path.join(Path, subject + '_flirt_defField.nii.gz')
    defField_magnitude_path = os.path.join(Path, subject + '_flirt_defField_magnitude.nii.gz')

    or_image = nib.load(or_path)

    with open(flirt_path, 'r') as f:
        txt = f.read()
        flirt = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
    print(flirt)
    with open(invflirt_path, 'r') as f:
        txt = f.read()
        invflirt = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
    print(invflirt)
    with open(flirtRigid_path, 'r') as f:
        txt = f.read()
        flirtRigid = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
    print(flirtRigid)
    with open(invflirtRigid_path, 'r') as f:
        txt = f.read()
        invflirtRigid = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
    print(invflirtRigid)

    # Calculate inverse flirt and check whether it matches with read inverse flirt

    #
    defField = imageFLIRT2defField(ref_image, invflirt)
    defField_image = nib.Nifti1Image(defField, ref_image.affine)
    nib.save(defField_image, defField_path)
    defField_magnitude = np.sqrt(np.sum(np.power(defField, 2), axis=3))
    defField_magnitude_image = nib.Nifti1Image(defField_magnitude, ref_image.affine)
    nib.save(defField_magnitude_image, defField_magnitude_path)

    # Assign syn defFIeld_magnitude's to 4D array
    defField_magnitude_flirt[:,:,:,iData] = defField_magnitude



    # Load SyN of subject and get inverse field
    defField_path = os.path.join(Path, subject + '_flirted_syn_defField.nii.gz')
    defField_magnitude_path = os.path.join(Path, subject + '_flirted_syn_defField_magnitude.nii.gz')
    syn_path = os.path.join(Path, subject + '_flirted_syn.pickle.gz')
    with open(syn_path, 'rb') as f:
        [mapping, metric, level_iters, sdr] = load(f, compression='gzip')

    defField = mapping.get_backward_field()
    del mapping, metric, level_iters, sdr

    defField_image = nib.Nifti1Image(defField, ref_image.affine)
    nib.save(defField_image, defField_path)
    defField_magnitude = np.sqrt(np.sum(np.power(defField, 2), axis=3))
    defField_magnitude_image = nib.Nifti1Image(defField_magnitude, ref_image.affine)
    nib.save(defField_magnitude_image, defField_magnitude_path)

    # Assign syn defFIeld_magnitude's to 4D array
    defField_magnitude_syn[:,:,:,iData] = defField_magnitude

# Average defField magnitude per genotype
genotype_WT_logical = np.array([genotype == 'WT' for genotype in genotype_list])
genotype_KO_logical = np.logical_not(genotype_WT_logical)
defField_magnitude_flirt_meanWT = np.empty(defField_magnitude.shape[0:3])
defField_magnitude_syn_meanWT = np.empty(defField_magnitude.shape[0:3])
defField_magnitude_flirt_meanKO = np.empty(defField_magnitude.shape[0:3])
defField_magnitude_syn_meanKO = np.empty(defField_magnitude.shape[0:3])
defField_magnitude_flirt_pval = np.empty(defField_magnitude.shape[0:3])
defField_magnitude_syn_pval = np.empty(defField_magnitude.shape[0:3])
for i in range(defField_magnitude.shape[0]):
    for j in range(defField_magnitude.shape[1]):
        for k in range(defField_magnitude.shape[2]):
            defField_magnitude_flirt_genotype = np.squeeze(defField_magnitude_flirt[i,j,k,:])
            defField_magnitude_syn_genotype = np.squeeze(defField_magnitude_syn[i,j,k,:])

            defField_magnitude_flirt_meanWT[i,j,k] = np.mean(defField_magnitude_flirt_genotype[genotype_WT_logical])
            defField_magnitude_syn_meanWT[i,j,k]   = np.mean(defField_magnitude_syn_genotype[  genotype_WT_logical])
            defField_magnitude_flirt_meanKO[i,j,k] = np.mean(defField_magnitude_flirt_genotype[genotype_KO_logical])
            defField_magnitude_syn_meanKO[i,j,k]   = np.mean(defField_magnitude_syn_genotype[  genotype_KO_logical])

            [_, defField_magnitude_flirt_pval[i,j,k]] = ttest_ind(defField_magnitude_flirt_genotype[genotype_WT_logical],
                                                                  defField_magnitude_flirt_genotype[genotype_KO_logical],
                                                                  equal_var=False)
            [_, defField_magnitude_syn_pval[i,j,k]] = ttest_ind(defField_magnitude_syn_genotype[genotype_WT_logical],
                                                                defField_magnitude_syn_genotype[genotype_KO_logical],
                                                                equal_var=False)
