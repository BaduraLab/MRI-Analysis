import numpy as np
from scipy import spatial
import numpy_indexed as npi
import nibabel as nib
import PIL.Image
import os
import glob
import pandas as pd


def imageAdjustNN(input, input_logical, correction, correction_logical):
    # Grid for point generation
    X, Y, Z = np.mgrid[0:input.shape[0]:1, 0:input.shape[1]:1, 0:input.shape[2]:1]

    # Points
    input_points = np.vstack((X[input_logical],
                              Y[input_logical],
                              Z[input_logical])).transpose()  # Get old automatic points
    correction_points = np.vstack((X[correction_logical],
                                   Y[correction_logical],
                                   Z[correction_logical])).transpose()  # Get new manual points

    # Correction tree
    correction_tree = spatial.KDTree(correction_points)  # Get old automatic tree

    # Go through input points and interpolate value to nearest correction value
    for iInput in range(input_points.shape[0]):
        add_index = tuple(input_points[iInput, :])
        closest_annotated_index = tuple(correction_points[correction_tree.query(add_index)[1]])
        input[add_index] = correction[closest_annotated_index]

    return input


def remap_3D(input, from_list, to_list):
    input = np.round(input)  # always annotation so should never be non-integer
    input = input.astype(int)  # always annotation so should never be non-integer

    input_shape = input.shape
    input = input.reshape(-1)
    input = npi.remap(input, from_list, to_list)
    input = input.reshape(input_shape)

    return input


def save_image(output, input_template_image, output_path):
    output_image = nib.Nifti1Image(output,
                                   input_template_image.affine,
                                   input_template_image.header)
    nib.save(output_image, output_path)


def reorient_image(input_path):
    input_image = nib.load(input_path)

    print(nib.aff2axcodes(input_image.affine))

    output_image = nib.as_closest_canonical(input_image)
    output_path = input_path.split('.')[0] + '_reoriented.nii.gz'
    print(output_path)

    print(nib.aff2axcodes(output_image.affine))

    nib.save(output_image, output_path)


def imageFolder2gif(folder_path, output_gif_path):
    filepath_list = glob.glob(os.path.join(folder_path, '*'))
    im_list = [PIL.Image.open(filepath) for filepath in filepath_list]
    im_list[0].save(output_gif_path,
                    ave_all=True,
                    append_images=im_list[1:],
                    optimize=False, duration=10, loop=0)


# Function to compute volumes for image
def subjectPath2volumeTable(subject_path):
    # Compute voxel numbers and volumes and output to table

    # Load image
    subject_image = nib.load(subject_path)
    subject = subject_image.get_fdata()

    # Get voxel volume
    print(subject_image.header['pixdim'][1:4])
    voxel_volume = np.prod(subject_image.header['pixdim'][1:4])  # should be in mm^3
    print('voxel volume = ' + str(voxel_volume) + 'mm^3')

    # Calculate volumes
    [iVoxel, nVoxel] = np.unique(np.int64(np.round(subject)),
                                 return_counts=True)
    vVoxel = nVoxel * voxel_volume
    print('total volume = ' + str(np.sum(vVoxel[iVoxel != 0])) + 'mm^3')

    # Output to DataFrame
    volume_table = pd.DataFrame(
        {'VolumeInteger': iVoxel,
         'VoxelNumber': nVoxel,
         'Volume': vVoxel})

    return volume_table


# Function to compute volumes for image
def imageFLIRT2defField(image, flirt):
    M_image = image.affine[:3, :3]
    abc_image = image.affine[:3, 3]
    M_flirt = flirt[:3, :3]
    abc_flirt = flirt[:3, 3]

    image_shape = image.shape  #######################
    defField = np.empty(shape=(image_shape[0], image_shape[1], image_shape[2], 3))
    for i in range(image_shape[0]):
        for j in range(image_shape[1]):
            for k in range(image_shape[2]):
                posVec = M_image.dot([i, j, k]) + abc_image
                flirtVec = M_flirt.dot(posVec) + abc_flirt
                defVec = flirtVec - posVec
                defField[i, j, k, :] = defVec

    return defField


# Function to zeropad 3D array
def zeroPadImage(input_3D_numpy, padRatio):
    nonzero_logical = np.where(input_3D_numpy > 0)
    edgeMax = np.max(nonzero_logical, axis=1)
    edgeMin = np.min(nonzero_logical, axis=1)
    edgeMax = edgeMax+1
    edgeDis = edgeMax - edgeMin
    edgeDis_max = np.max(edgeDis)
    padLength = int(np.round(edgeDis_max * padRatio))

    print(edgeMax)
    print(edgeMin)
    print(edgeDis_max)
    print(padLength)

    output = input_3D_numpy[edgeMin[0]:edgeMax[0], edgeMin[1]:edgeMax[1], edgeMin[2]:edgeMax[2]]

    # check whether amount of nonzero elements is still equal after cropping
    print((f'Original number of nonzero elemenets = {np.sum(input_3D_numpy>0)}'))
    print((f'Maximally cropped output number of nonzero elemenets = {np.sum(output>0)}'))

    output = np.pad(output,
                    ((padLength, padLength),
                     (padLength, padLength),
                     (padLength, padLength)),
                    'constant')
    print((f'Zero padded output number of nonzero elemenets = {np.sum(output>0)}'))

    return output
