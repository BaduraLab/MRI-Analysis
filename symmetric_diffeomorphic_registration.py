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



hardi_fname, hardi_bval_fname, hardi_bvec_fname = get_fnames('stanford_hardi')

stanford_b0, stanford_b0_affine = load_nifti(hardi_fname)
stanford_b0 = np.squeeze(stanford_b0)[..., 0]



t1_fname, b0_fname = get_fnames('syn_data')
syn_b0, syn_b0_affine = load_nifti(b0_fname)



from dipy.segment.mask import median_otsu
stanford_b0_masked, stanford_b0_mask = median_otsu(stanford_b0,
                                                   median_radius=4,
                                                   numpass=4)
syn_b0_masked, syn_b0_mask = median_otsu(syn_b0, median_radius=4, numpass=4)

static = stanford_b0_masked
static_affine = stanford_b0_affine
moving = syn_b0_masked
moving_affine = syn_b0_affine



# FLIRT
pre_align = np.array([[1.02783543e+00, -4.83019053e-02, -6.07735639e-02, -2.57654118e+00],
                      [4.34051706e-03, 9.41918267e-01, -2.66525861e-01, 3.23579799e+01],
                      [5.34288908e-02, 2.90262026e-01, 9.80820307e-01, -1.46216651e+01],
                      [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00]])



from dipy.align.imaffine import AffineMap
affine_map = AffineMap(pre_align,
                       static.shape, static_affine,
                       moving.shape, moving_affine)

resampled = affine_map.transform(moving)



regtools.overlay_slices(static, resampled, None, 1, 'Static', 'Moving',
                        'input_3d.png')



metric = CCMetric(3)



level_iters = [10, 10, 5]
sdr = SymmetricDiffeomorphicRegistration(metric, level_iters)



mapping = sdr.optimize(static, moving, static_affine, moving_affine, pre_align)



warped_moving = mapping.transform(moving)



regtools.overlay_slices(static, warped_moving, None, 1, 'Static',
                        'Warped moving', 'warped_moving.png')



warped_static = mapping.transform_inverse(static)
regtools.overlay_slices(warped_static, moving, None, 1, 'Warped static',
                        'Moving', 'warped_static.png')