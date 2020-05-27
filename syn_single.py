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
import plotly.graph_objects as go



# Define
reference_path = '/usr/local/fsl/data/standard/allen_new'
input_path = '/home/enzo/Desktop/Data/Mouse/Processed_New/WT_50'
FLIRT_folder_path = os.path.join(input_path, 'FLIRT')
FLIRT_path = glob.glob(FLIRT_folder_path+'/*warpaffine.mat')
reference_image_path = os.path.join(reference_path, 'average_template_25_to_AMBMC_flirted.nii.gz')
annotation_image_path = os.path.join(reference_path, 'annotation_25_to_AMBMC_flirted.nii.gz')
input_image_path = os.path.join(FLIRT_folder_path, 'WT_50_flirted.nii.gz')



# Load images
reference_image = nib.load(reference_image_path)
reference = reference_image.get_fdata()
annotation_image = nib.load(annotation_image_path)
annotation = annotation_image.get_fdata()
input_image = nib.load(input_image_path)
input = input_image.get_fdata()
# FLIRT = scipy.io.loadmat(FLIRT_path[0])
FLIRT = np.array([[0.8719689183,  -0.003598787672,  0.01898341971,  -1.128021769],
                        [0.00478144699,  0.8647957243,  0.1031076254,  3.502086638  ],
                        [0.01072608439,  -0.09020316308,  0.9738558691,  0.677238873  ],
                    [0,  0,  0,  1  ]])


# syn
metric = CCMetric(3)
level_iters = [10, 10, 5]
sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)

mapping = sdr.optimize(static=reference, moving=input,
                       static_grid2world=reference_image.get_qform(), moving_grid2world=input_image.get_qform())

warped_image = mapping.transform(input)

warped_reference = mapping.transform_inverse(reference)

warped_annotation = mapping.transform_inverse(annotation)

# evaluate warps
X, Y, Z = np.mgrid[0:input.shape[0]:1,650:651:1,0:input.shape[2]:1]

fig1 = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=input[0:input.shape[0]:1,650:651:1,0:input.shape[2]:1].flatten(),
    opacity=1, # needs to be small to see through all surfaces
    surface_count=17, # needs to be a large number for good volume rendering
    ))
fig1.show()

X, Y, Z = np.mgrid[0:reference.shape[0]:1,650:651:1,0:reference.shape[2]:1]

fig2 = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=reference[0:reference.shape[0]:1,650:651:1,0:reference.shape[2]:1].flatten(),
    opacity=1, # needs to be small to see through all surfaces
    surface_count=17, # needs to be a large number for good volume rendering
    ))
fig2.show()

X, Y, Z = np.mgrid[0:warped_image.shape[0]:1,650:651:1,0:warped_image.shape[2]:1]

fig3 = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=warped_image[0:warped_image.shape[0]:1,650:651:1,0:warped_image.shape[2]:1].flatten(),
    opacity=1, # needs to be small to see through all surfaces
    surface_count=17, # needs to be a large number for good volume rendering
    ))
fig3.show()

X, Y, Z = np.mgrid[0:warped_annotation.shape[0]:1,650:651:1,0:warped_annotation.shape[2]:1]

fig4 = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=warped_annotation[0:warped_annotation.shape[0]:1,650:651:1,0:warped_annotation.shape[2]:1].flatten(),
    opacity=1, # needs to be small to see through all surfaces
    surface_count=17, # needs to be a large number for good volume rendering
    ))
fig4.show()
