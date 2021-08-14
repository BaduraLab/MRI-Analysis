# Get native coordinates of human and mouse patient and WT images
import os
import nibabel as nib

# Define
import numpy as np
import pandas as pd

from functions import fslCoords2voxCoords, voxCoords2fslCoords
from scipy.ndimage.measurements import center_of_mass
from functions import save_image
import numpy as np
from scipy import signal
import pandas as pd





# # first build the smoothing kernel
# sigma = 1.0     # width of kernel
# x = np.arange(-3,4,1)   # coordinate arrays -- make sure they contain 0!
# y = np.arange(-3,4,1)
# z = np.arange(-3,4,1)
# xx, yy, zz = np.meshgrid(x,y,z)
# kernel = np.exp(-(xx**2 + yy**2 + zz**2)/(2*sigma**2))
#
# # apply to sample data
# data = np.zeros((11,11,11))
# data[5,5,5] = 5.
# filtered = signal.convolve(data, kernel, mode="same")
#
# # # check output
# # print filtered[:,5,5]
#
# mouse Coords = np.array([4.96, -8.38, 3.82]) # SN_CE [4.95, -8.375, 3.825]
# human Coords = np.array([10, -16, -13])  # still in flirtedRigid, SN_CE




# Define paths
analysis_path = os.path.join('Data', 'Analysis')
data_human_path = os.path.join('Data', 'Human', 'Processed')
data_mouse_path = os.path.join('Data', 'Mouse', 'Processed_Old')

# Mouse control
mouse_control_subject_path = os.path.join(data_mouse_path, 'WT_50')
mouse_control_native = os.path.join(data_mouse_path, 'WT_50', 'FLIRT', 'WT_50_inmasked.nii.gz')
mouse_control_flirtedRigid = os.path.join(data_mouse_path, 'WT_50', 'retransform', 'WT_50_flirtedRigid.nii.gz')
mouse_control_flirted = os.path.join(data_mouse_path, 'WT_50', 'retransform', 'WT_50_flirted.nii.gz')
mouse_control_flirtRigid = os.path.join(data_mouse_path, 'WT_50', 'retransform', 'flirtRigid.mat')
mouse_control_flirt = os.path.join(data_mouse_path, 'WT_50', 'retransform', 'flirt.mat')
mouse_control_invflirtRigid = os.path.join(data_mouse_path, 'WT_50', 'retransform', 'invflirtRigid.mat')
mouse_control_invflirt = os.path.join(data_mouse_path, 'WT_50', 'retransform', 'invflirt.mat')
mouse_control_ref_pointSignal = os.path.join(data_mouse_path, 'WT_50', 'retransform', 'WT_50_pointSignal_ref_SN_CE.nii.gz')
mouse_control_flirtedRigid_pointSignal = os.path.join(data_mouse_path, 'WT_50', 'retransform', 'WT_50_pointSignal_flirtedRigid_SN_CE.nii.gz')
mouse_control_native_pointSignal = os.path.join(data_mouse_path, 'WT_50', 'retransform', 'WT_50_pointSignal_native_SN_CE.nii.gz')
mouse_control_Coords = np.array([4.96, -8.38, 3.82]) # SN_CE
os.system('convert_xfm -omat ' + mouse_control_invflirtRigid + ' -inverse ' + mouse_control_flirtRigid)
with open(mouse_control_invflirtRigid, 'r') as f:
    txt = f.read()
    mouse_control_invflirtRigid_mat = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
print(mouse_control_invflirtRigid_mat)
os.system('convert_xfm -omat ' + mouse_control_invflirt + ' -inverse ' + mouse_control_flirt)
with open(mouse_control_invflirt, 'r') as f:
    txt = f.read()
    mouse_control_invflirt_mat = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
print(mouse_control_invflirt_mat)

# Mouse KO
mouse_KO_subject_path = os.path.join(data_mouse_path, 'KO_70')
mouse_KO_native = os.path.join(data_mouse_path, 'KO_70', 'FLIRT', 'KO_70_inmasked.nii.gz')
mouse_KO_flirtedRigid = os.path.join(data_mouse_path, 'KO_70', 'retransform', 'KO_70_flirtedRigid.nii.gz')
mouse_KO_flirted = os.path.join(data_mouse_path, 'KO_70', 'retransform', 'KO_70_flirted.nii.gz')
mouse_KO_flirtRigid = os.path.join(data_mouse_path, 'KO_70', 'retransform', 'flirtRigid.mat')
mouse_KO_flirt = os.path.join(data_mouse_path, 'KO_70', 'retransform', 'flirt.mat')
mouse_KO_invflirtRigid = os.path.join(data_mouse_path, 'KO_70', 'retransform', 'invflirtRigid.mat')
mouse_KO_invflirt = os.path.join(data_mouse_path, 'KO_70', 'retransform', 'invflirt.mat')
mouse_KO_ref_pointSignal = os.path.join(data_mouse_path, 'KO_70', 'retransform', 'KO_70_pointSignal_ref_SN_CE.nii.gz')
mouse_KO_flirtedRigid_pointSignal = os.path.join(data_mouse_path, 'KO_70', 'retransform', 'KO_70_pointSignal_flirtedRigid_SN_CE.nii.gz')
mouse_KO_native_pointSignal = os.path.join(data_mouse_path, 'KO_70', 'retransform', 'KO_70_pointSignal_native_SN_CE.nii.gz')
mouse_KO_Coords = np.array([4.96, -8.38, 3.82]) # SN_CE
os.system('convert_xfm -omat ' + mouse_KO_invflirtRigid + ' -inverse ' + mouse_KO_flirtRigid)
with open(mouse_KO_invflirtRigid, 'r') as f:
    txt = f.read()
    mouse_KO_invflirtRigid_mat = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
print(mouse_KO_invflirtRigid_mat)
os.system('convert_xfm -omat ' + mouse_KO_invflirt + ' -inverse ' + mouse_KO_flirt)
with open(mouse_KO_invflirt, 'r') as f:
    txt = f.read()
    mouse_KO_invflirt_mat = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
print(mouse_KO_invflirt_mat)

# Human control
human_control_subject_path = os.path.join(data_human_path, 'control4')
human_control_native = os.path.join(data_human_path, 'control4', 'control4_reoriented.nii.gz')
human_control_flirtedRigid = os.path.join(data_human_path, 'control4', 'retransform', 'control4_flirtedRigid.nii.gz')
human_control_flirted = os.path.join(data_human_path, 'control4', 'retransform', 'control4_flirted.nii.gz')
human_control_flirtRigid = os.path.join(data_human_path, 'control4', 'retransform', 'flirtRigid.mat')
human_control_flirt = os.path.join(data_human_path, 'control4', 'retransform', 'flirt.mat')
human_control_invflirtRigid = os.path.join(data_human_path, 'control4', 'retransform', 'invflirtRigid.mat')
human_control_invflirt = os.path.join(data_human_path, 'control4', 'retransform', 'invflirt.mat')
human_control_ref_pointSignal = os.path.join(data_human_path, 'control4', 'retransform', 'control4_pointSignal_ref_SN_CE.nii.gz')
human_control_flirtedRigid_pointSignal = os.path.join(data_human_path, 'control4', 'retransform', 'control4_pointSignal_flirtedRigid_SN_CE.nii.gz')
human_control_native_pointSignal = os.path.join(data_human_path, 'control4', 'retransform', 'control4_pointSignal_native_SN_CE.nii.gz')
human_control_Coords = np.array([10, -16, -13])  # still in flirtedRigid, SN_CE
os.system('convert_xfm -omat ' + human_control_invflirtRigid + ' -inverse ' + human_control_flirtRigid)
with open(human_control_invflirtRigid, 'r') as f:
    txt = f.read()
    human_control_invflirtRigid_mat = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
print(human_control_invflirtRigid_mat)
os.system('convert_xfm -omat ' + human_control_invflirt + ' -inverse ' + human_control_flirt)
with open(human_control_invflirt, 'r') as f:
    txt = f.read()
    human_control_invflirt_mat = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
print(human_control_invflirt_mat)

# Human patient
human_patient_subject_path = os.path.join(data_human_path, 'patient')
human_patient_native = os.path.join(data_human_path, 'patient', 'patient_reoriented.nii.gz')
human_patient_flirtedRigid = os.path.join(data_human_path, 'patient', 'retransform', 'patient_flirtedRigid.nii.gz')
human_patient_flirted = os.path.join(data_human_path, 'patient', 'retransform', 'patient_flirted.nii.gz')
human_patient_flirtRigid = os.path.join(data_human_path, 'patient', 'retransform', 'flirtRigid.mat')
human_patient_flirt = os.path.join(data_human_path, 'patient', 'retransform', 'flirt.mat')
human_patient_invflirtRigid = os.path.join(data_human_path, 'patient', 'retransform', 'invflirtRigid.mat')
human_patient_invflirt = os.path.join(data_human_path, 'patient', 'retransform', 'invflirt.mat')
human_patient_ref_pointSignal = os.path.join(data_human_path, 'patient', 'retransform', 'patient_pointSignal_ref_SN_CE.nii.gz')
human_patient_flirtedRigid_pointSignal = os.path.join(data_human_path, 'patient', 'retransform', 'patient_pointSignal_flirtedRigid_SN_CE.nii.gz')
human_patient_native_pointSignal = os.path.join(data_human_path, 'patient', 'retransform', 'patient_pointSignal_native_SN_CE.nii.gz')
human_patient_Coords = np.array([10, -16, -13])  # still in flirtedRigid, SN_CE
os.system('convert_xfm -omat ' + human_patient_invflirtRigid + ' -inverse ' + human_patient_flirtRigid)
with open(human_patient_invflirtRigid, 'r') as f:
    txt = f.read()
    human_patient_invflirtRigid_mat = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
print(human_patient_invflirtRigid_mat)
os.system('convert_xfm -omat ' + human_patient_invflirt + ' -inverse ' + human_patient_flirt)
with open(human_patient_invflirt, 'r') as f:
    txt = f.read()
    human_patient_invflirt_mat = np.array([[float(num) for num in item.split()] for item in txt.split('\n')[:-1]])
print(human_patient_invflirt_mat)


# Load image
# Get coords flirted by using qform and defined voxel coordinates
# Convert coords to native coords
# Convert native coords to voxel coords
## function

# def ref_voxCoords_convert(ref_voxCoords, ref_image, invflirt, flirtedRigid_image, invflirtRigid, native_image):
#     ref_Q = ref_image.get_qform()
#     dims = np.abs(np.array([ref_Q[0, 0], ref_Q[1, 1], ref_Q[2, 2]]))
#
#     # Get dims_native by assuming uniform spatial sampling in native space - not true in case of human images...
#     dims_native = np.power(np.abs(np.linalg.det(ref_image.get_qform())), 1/3)
#     dims_native = np.array([dims_native, dims_native, dims_native])
#
#     ref_Coords = np.matmul(ref_Q, ref_voxCoords)
#
#     ref_fslCoords = voxCoords2fslCoords(ref_voxCoords, ref_image, dims)
#     flirtedRigid_fslCoords = np.matmul(invflirt, ref_fslCoords) ###
#     flirtedRigid_voxCoords = fslCoords2voxCoords(flirtedRigid_fslCoords, ref_image, dims)
#     # flirtedRigid_invQ = np.linalg.inv(flirtedRigid_image.get_qform())
#     flirtedRigid_Coords = np.matmul(flirtedRigid_image.get_qform(), flirtedRigid_voxCoords)
#
#     # flirtedRigid_fslCoords = voxCoords2fslCoords(flirtedRigid_voxCoords, ref_image, dims_native)
#     native_fslCoords = np.matmul(invflirtRigid, flirtedRigid_fslCoords) ###
#     native_voxCoords = fslCoords2voxCoords(native_fslCoords, ref_image, dims) #?# does this make sense? native to native fsl to flirtedRigid
#     # native_invQ = np.linalg.inv(native_image.get_qform())
#     native_Coords = np.matmul(native_image.get_qform(), native_voxCoords)
#
#     return (ref_Coords, flirtedRigid_Coords, native_Coords), \
#            (ref_voxCoords, flirtedRigid_voxCoords, native_voxCoords), \
#            (ref_fslCoords, flirtedRigid_fslCoords, native_fslCoords)
# need to convert ref_voxCoords to
def ref_voxCoords_convert(ref_Coords, ref_image, ref_pointSignal_path,
                          flirtedRigid_path, flirtedRigid_pointSignal_path, invflirtRigid_path,
                          native_path, native_pointSignal_path, invflirt_path):
    # Define point signal image
    ref_invQ = np.linalg.inv(ref_image.get_qform())
    ref_voxCoords = np.matmul(ref_invQ, np.concatenate([ref_Coords, [1]]))
    ref_voxCoords = ref_voxCoords[0:3].astype(int)
    ref_pointSignal = np.zeros(ref_image.shape)
    # ref_voxCoords[0] = ref_image.shape[0] - 1 - ref_voxCoords[0] # adjust ref_voxCoords, unknown whether this is always necessary, maybe only if np.linalg.det(image.get_qform())>0?
    ref_pointSignal[tuple(ref_voxCoords)] = 1
    # for i in range(3):
    # ref_pointSignal = signal.convolve(ref_pointSignal, kernel, mode="same")
    ref_COM_voxCoords = center_of_mass(ref_pointSignal)
    print(f'ref_COM_voxCoords = {ref_COM_voxCoords}')
    save_image(ref_pointSignal, ref_image, ref_pointSignal_path)

    # # Get reference Coords
    # ref_Coords = np.matmul(ref_image.get_qform(), np.concatenate([ref_voxCoords, [1]]))

    # Transform point signal image to flirtedRigid space, calculate COM
    print('flirt apply invflirt start')
    os.system('flirt -in ' + ref_pointSignal_path + ' \
                     -ref ' + flirtedRigid_path + ' \
                     -out ' + flirtedRigid_pointSignal_path + ' \
                     -init ' + invflirt_path + ' \
                     -applyxfm \
                     -verbose 0')
    flirtedRigid_pointSignal_image = nib.load(flirtedRigid_pointSignal_path)
    flirtedRigid_pointSignal = flirtedRigid_pointSignal_image.get_fdata()
    flirtedRigid_COM_voxCoords = center_of_mass(flirtedRigid_pointSignal)
    flirtedRigid_COM_Coords = np.matmul(flirtedRigid_pointSignal_image.get_qform(), np.concatenate([flirtedRigid_COM_voxCoords, [1]]))
    flirtedRigid_COM_Coords = flirtedRigid_COM_Coords[0:3]


    # Transform point signal image to native space, calculate COM
    print('flirt apply invflirtRigid start')
    os.system('flirt -in ' + flirtedRigid_pointSignal_path + ' \
                     -ref ' + native_path + ' \
                     -out ' + native_pointSignal_path + ' \
                     -init ' + invflirtRigid_path + ' \
                     -applyxfm \
                     -verbose 0')
    native_pointSignal_image = nib.load(native_pointSignal_path)
    native_pointSignal = native_pointSignal_image.get_fdata()
    native_COM_voxCoords = center_of_mass(native_pointSignal)
    native_COM_Coords = np.matmul(native_pointSignal_image.get_qform(), np.concatenate([native_COM_voxCoords, [1]]))
    native_COM_Coords = native_COM_Coords[0:3]

    return (ref_voxCoords, flirtedRigid_COM_voxCoords, native_COM_voxCoords), \
           (ref_Coords, flirtedRigid_COM_Coords, native_COM_Coords)

# Human patient
(human_patient_voxCoords_tuple, human_patient_Coords_tuple) = ref_voxCoords_convert(ref_Coords=human_patient_Coords,
                                                                        ref_image=nib.load(human_patient_flirted),
                                                                        ref_pointSignal_path=human_patient_ref_pointSignal,
                                                                        flirtedRigid_path=human_patient_flirtedRigid,
                                                                        flirtedRigid_pointSignal_path=human_patient_flirtedRigid_pointSignal,
                                                                        invflirtRigid_path=human_patient_invflirtRigid,
                                                                        native_path=human_patient_native,
                                                                        native_pointSignal_path=human_patient_native_pointSignal,
                                                                        invflirt_path=human_patient_invflirt)

# Human control
(human_control_voxCoords_tuple, human_control_Coords_tuple) = ref_voxCoords_convert(ref_Coords=human_control_Coords,
                                                                        ref_image=nib.load(human_control_flirted),
                                                                        ref_pointSignal_path=human_control_ref_pointSignal,
                                                                        flirtedRigid_path=human_control_flirtedRigid,
                                                                        flirtedRigid_pointSignal_path=human_control_flirtedRigid_pointSignal,
                                                                        invflirtRigid_path=human_control_invflirtRigid,
                                                                        native_path=human_control_native,
                                                                        native_pointSignal_path=human_control_native_pointSignal,
                                                                        invflirt_path=human_control_invflirt)

# Mouse KO
(mouse_KO_voxCoords_tuple, mouse_KO_Coords_tuple) = ref_voxCoords_convert(ref_Coords=mouse_KO_Coords,
                                                                        ref_image=nib.load(mouse_KO_flirted),
                                                                        ref_pointSignal_path=mouse_KO_ref_pointSignal,
                                                                        flirtedRigid_path=mouse_KO_flirtedRigid,
                                                                        flirtedRigid_pointSignal_path=mouse_KO_flirtedRigid_pointSignal,
                                                                        invflirtRigid_path=mouse_KO_invflirtRigid,
                                                                        native_path=mouse_KO_native,
                                                                        native_pointSignal_path=mouse_KO_native_pointSignal,
                                                                        invflirt_path=mouse_KO_invflirt)

# Mouse control
(mouse_control_voxCoords_tuple, mouse_control_Coords_tuple) = ref_voxCoords_convert(ref_Coords=mouse_control_Coords,
                                                                        ref_image=nib.load(mouse_control_flirted),
                                                                        ref_pointSignal_path=mouse_control_ref_pointSignal,
                                                                        flirtedRigid_path=mouse_control_flirtedRigid,
                                                                        flirtedRigid_pointSignal_path=mouse_control_flirtedRigid_pointSignal,
                                                                        invflirtRigid_path=mouse_control_invflirtRigid,
                                                                        native_path=mouse_control_native,
                                                                        native_pointSignal_path=mouse_control_native_pointSignal,
                                                                        invflirt_path=mouse_control_invflirt)

# Create table
strID = ['human_patient_SN_CE', 'human_control_SN_CE', 'mouse_KO_SN_CE', 'mouse_control_SN_CE']
data_subject_path_list = [human_patient_subject_path, human_control_subject_path, mouse_KO_subject_path, mouse_control_subject_path]
# synned_path_list = [human_patient_synned[0], human_control_voxCoords_tuple[0], mouse_KO_voxCoords_tuple[0], mouse_control_voxCoords_tuple[0]]
# flirted_path_list
ref_voxCoords = [list(human_patient_voxCoords_tuple[0]), list(human_control_voxCoords_tuple[0]), list(mouse_KO_voxCoords_tuple[0]), list(mouse_control_voxCoords_tuple[0])]
flirtedRigid_voxCoords = [list(human_patient_voxCoords_tuple[1]), list(human_control_voxCoords_tuple[1]), list(mouse_KO_voxCoords_tuple[1]), list(mouse_control_voxCoords_tuple[1])]
native_voxCoords = [list(human_patient_voxCoords_tuple[2]), list(human_control_voxCoords_tuple[2]), list(mouse_KO_voxCoords_tuple[2]), list(mouse_control_voxCoords_tuple[2])]
ref_Coords = [list(human_patient_Coords_tuple[0]), list(human_control_Coords_tuple[0]), list(mouse_KO_Coords_tuple[0]), list(mouse_control_Coords_tuple[0])]
flirtedRigid_Coords = [list(human_patient_Coords_tuple[1]), list(human_control_Coords_tuple[1]), list(mouse_KO_Coords_tuple[1]), list(mouse_control_Coords_tuple[1])]
native_Coords = [list(human_patient_Coords_tuple[2]), list(human_control_Coords_tuple[2]), list(mouse_KO_Coords_tuple[2]), list(mouse_control_Coords_tuple[2])]
coords_table = pd.DataFrame({'strID': strID,
                             'data_subject_path': data_subject_path_list,
                             'ref_voxCoords': ref_voxCoords,
                             'flirtedRigid_voxCoords': flirtedRigid_voxCoords,
                             'native_voxCoords': native_voxCoords,
                             'ref_Coords': ref_Coords,
                             'flirtedRigid_Coords': flirtedRigid_Coords,
                             'native_Coords': native_Coords})
coords_table.to_csv(os.path.join(analysis_path, 'coords_table.csv'))
# data_subject_path




# # Human control
# human_control_voxCoords = np.array([201, 218, 124, 1])  # still in flirtedRigid
# (human_control_Coords, human_control_voxCoords, human_control_fslCoords) = ref_voxCoords_convert(ref_voxCoords=human_control_voxCoords,
#                                                                         ref_image=human_control_flirted,
#                                                                         invflirt=human_control_invflirt_mat,
#                                                                         flirtedRigid_image=human_control_flirtedRigid,
#                                                                         invflirtRigid=human_control_invflirtRigid_mat,
#                                                                         native_image=human_control_native)
#
# # Mouse KO
# mouse_KO_voxCoords = np.array([286, 442, 89, 1])
# (mouse_KO_Coords, mouse_KO_voxCoords, mouse_KO_fslCoords) = ref_voxCoords_convert(ref_voxCoords=mouse_KO_voxCoords,
#                                                               ref_image=mouse_KO_flirted,
#                                                               invflirt=mouse_KO_invflirt_mat,
#                                                               flirtedRigid_image=mouse_KO_flirtedRigid,
#                                                               invflirtRigid=mouse_KO_invflirtRigid_mat,
#                                                               native_image=mouse_KO_native)
#
# # Mouse WT
# mouse_control_voxCoords = np.array([286, 442, 89, 1])
# (mouse_control_Coords, mouse_control_voxCoords, mouse_control_fslCoords) = ref_voxCoords_convert(ref_voxCoords=mouse_control_voxCoords,
#                                                                         ref_image=mouse_control_flirted,
#                                                                         invflirt=mouse_control_invflirt_mat,
#                                                                         flirtedRigid_image=mouse_control_flirtedRigid,
#                                                                         invflirtRigid=mouse_control_invflirtRigid_mat,
#                                                                         native_image=mouse_control_native)
#
# # mricrogl
# # make images with mricrogl at points (have both vox coords and coords)
# # multi background
# # multi full overlay CE and SN (CE not that necessary for mouse)
