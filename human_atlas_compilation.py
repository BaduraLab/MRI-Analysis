import os
import nibabel as nib
import datetime
import numpy as np
from dipy.align.imwarp import SymmetricDiffeomorphicRegistration
from dipy.align.imwarp import DiffeomorphicMap
from dipy.align.metrics import CCMetric

# Define paths
fsl_path = '/usr/local/fsl'
data_path = '/home/enzo/Desktop/Data'
# annotation_path = atlas_fsl_dir+'/data/atlases/Talairach/Talairach-labels-1mm.nii.gz'
annotation_path = fsl_path+'/data/atlases/Cerebellum/Cerebellum-MNIfnirt-maxprob-thr50-1mm.nii.gz'
template_path = fsl_path+'/data/standard/MNI152_T1_1mm_brain.nii.gz'
# input_path = data_path+'/Human/controls/MPRAGE_average_deface_contr1.nii' # input should be only brain, no skull
# input_path = data_path+'/Human/controls/MPRAGE_average_deface_contr2.nii' # input should be only brain, no skull
input_path = data_path+'/Human/Pax5_mutation/MPRAGE_average_deface.nii' # input should be only brain, no skull
input_flirted_path = input_path.split('.')[0]+'_flirted.nii.gz'
input_flirt_path = input_path.split('.')[0]+'_flirt.mat'
input_invflirt_path = input_path.split('.')[0]+'_invflirt.mat'
input_flirted_synned_path = input_path.split('.')[0]+'_flirted_synned.nii.gz'
annotation_invsynned_path = input_path.split('.')[0]+'_flirted_annotation.nii.gz'
annotation_invsynned_invflirted_path = input_path.split('.')[0]+'_annotation.nii.gz'



# FLIRT
# Invert FLIRT warped annotation back to subject space
os.system('flirt -in ' + input_path + ' \
                 -ref ' + template_path + ' \
                 -out ' + input_flirted_path + ' \
                 -omat ' + input_flirt_path + ' \
                 -verbose 1')




# Load images
annotation_image = nib.load(annotation_path)
annotation = annotation_image.get_fdata()
template_image = nib.load(template_path)
template = template_image.get_fdata()
input_flirted_image = nib.load(input_flirted_path)
input_flirted = input_flirted_image.get_fdata()



# syn
metric = CCMetric(3)
level_iters = [10, 10, 5, 5, 5]
print(datetime.datetime.now())
sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)

mapping = sdr.optimize(static=template, moving=input_flirted,
                       static_grid2world=template_image.get_qform(), moving_grid2world=input_flirted_image.get_qform())

input_flirted_synned = mapping.transform(input_flirted)

annotation_invsynned = mapping.transform_inverse(annotation, interpolation='nearest')

# save some (inv)synned warped images, mostly to check if everything went well
input_flirted_synned_image = nib.Nifti1Image(input_flirted_synned, np.eye(4))
input_flirted_synned_image.set_qform(template_image.get_qform(), code=1)
input_flirted_synned_image.set_sform(np.eye(4), code=0)
nib.save(input_flirted_synned_image, input_flirted_synned_path)

annotation_invsynned_image = nib.Nifti1Image(annotation_invsynned, np.eye(4))
annotation_invsynned_image.set_qform(input_flirted_image.get_qform(), code=1)
annotation_invsynned_image.set_sform(np.eye(4), code=0)
nib.save(annotation_invsynned_image, annotation_invsynned_path)



# inflirt invsynned annotation to flirted image to get annotation of original image
os.system('convert_xfm -omat '+input_invflirt_path+' -inverse '+input_flirt_path)
os.system('flirt -in ' + annotation_invsynned_path + ' \
                 -ref ' + input_path + ' \
                 -out ' + annotation_invsynned_invflirted_path + ' \
                 -init ' + input_invflirt_path + ' \
                 -applyxfm \
                 -interp nearestneighbour \
                 -verbose 1')
