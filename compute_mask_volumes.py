import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind



# Define
data_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse/Processed_New'
# voxel_volume = pow(0.05, 3)
voxel_volume = 0.000125

mouse_list = os.listdir(data_path)
nMouse = len(mouse_list)



mouse_voxel_number = mouse_list.copy()
mouse_mask_volume = mouse_list.copy()
for iMouse, Mouse in enumerate(mouse_list):
    mouse_mask_image = nib.load(os.path.join(data_path, Mouse, (Mouse+'_mask_t=500_v=380_k=6.mask.nii.gz')))
    mouse_mask_image_array = mouse_mask_image.get_fdata()

    mouse_voxel_number[iMouse] = np.sum(mouse_mask_image_array>0)
    mouse_mask_volume[iMouse] = mouse_voxel_number[iMouse]*voxel_volume

mouse_table = pd.DataFrame({'Mouse': mouse_list, 'MaskVoxelNumber': mouse_voxel_number, 'MaskVolume': mouse_mask_volume})
mouse_table['Genotype'] = mouse_table['Mouse'].str.split("_", n = 1, expand = True).iloc[:, 0]
mouse_table['Sex'] = mouse_table['Mouse'].str.split("_", n = 3, expand = True).iloc[:, 2]
mouse_table.loc[mouse_table['Sex']!='female', 'Sex'] = 'male'

plt.figure();
mouse_table[['MaskVolume', 'Genotype', 'Sex']].boxplot(by=['Genotype', 'Sex'])
plt.savefig('C:/Users/enzo/Downloads/Study/Current Courses/MEP/Pictures/Boxplot_MaskVolumes_ByGenotype')



mouse_table[mouse_table['Genotype']=='WT']['MaskVolume']
cat2 = mouse_table[mouse_table['Genotype']=='KO']

ttest_ind(mouse_table[mouse_table['Genotype']=='WT']['MaskVolume'],
          mouse_table[mouse_table['Genotype']=='KO']['MaskVolume'])
