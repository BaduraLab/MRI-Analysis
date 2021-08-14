import numpy as np
from scipy import spatial
import numpy_indexed as npi
import nibabel as nib
import PIL.Image
import os
import glob
import pandas as pd
import fractions, math
from numba import jit
import time
from scipy import stats


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
    print(f'Saving {output_path}')
    output_image = nib.Nifti1Image(output,
                                   input_template_image.affine,
                                   input_template_image.header)
    nib.save(output_image, output_path)
    return output_image


def reorient_image(input_path):
    input_image = nib.load(input_path)

    print(nib.aff2axcodes(input_image.affine))

    output_image = nib.as_closest_canonical(input_image)
    output_path = input_path.split('.')[0] + '_reoriented.nii.gz'
    print(output_path)

    print(nib.aff2axcodes(output_image.affine))

    nib.save(output_image, output_path)


def imageFolder2gif(folder_path, output_gif_path=None, fps=10, max_gif_size=50):
    import imageio
    from pygifsicle import optimize
    from skimage.transform import resize
    import sys

    # List image paths in image sequence folder
    filepath_list = glob.glob(os.path.join(folder_path, '*'))

    # Read sorted images
    images = []
    ymin = []; xmin = []
    ymax = []; xmax = []
    sorted_degrees_indices = np.argsort([int(filepath.split('.')[-2].split('_')[-1]) for filepath in filepath_list])
    for iFilename in sorted_degrees_indices:
        filename = filepath_list[iFilename]
        print(filename)
        im = imageio.imread(filename)
        if np.sum(im) != 0: # only add im if it contains nonzeros
            images.append(im)
            nonzero_yxCoord = np.where(np.max(im, 2))

            # Add minimum and maximum y and x coordinates for cropping
            ymin.append(np.min(nonzero_yxCoord[0]))
            ymax.append(np.max(nonzero_yxCoord[0]))
            xmin.append(np.min(nonzero_yxCoord[1]))
            xmax.append(np.max(nonzero_yxCoord[1]))

    # Calculate outermost cropping box corners
    ymin = np.min(ymin)
    ymax = np.max(ymax) + 1
    xmin = np.min(xmin)
    xmax = np.max(xmax) + 1

    # Crop images
    images_cropped = []
    for im in images:
        images_cropped.append(im[ymin:ymax, xmin:xmax])

    # If cropped image object is larger than max_gif_size, resize it
    # (eventual gif will be much smaller than max_gif_size)
    images_cropped_size = np.sum([imC.nbytes for imC in images_cropped]) / 1e6
    if images_cropped_size > max_gif_size:
        print(f'Size of cropped images = {images_cropped_size} MB')
        # According to the size of images_cropped,
        # downsample the images to a limit of 10 MB so that the gif will be < 10 MB
        resize_factor = np.sqrt(max_gif_size / images_cropped_size)
        print(f'Resize factor = {resize_factor}')
        original_shape = images_cropped[0].shape # Cropped image should all have same shape
        resize_shape = [int(Dim) for Dim in np.array(original_shape) * resize_factor]
        resize_shape[2] = original_shape[2]
        print(f'Resize shape = {resize_shape}')
        for iIMC in range(len(images_cropped)):
            images_cropped[iIMC] = resize(images_cropped[iIMC], resize_shape, order=3)  # Bi-cubic interpolation
            # images_cropped[iIMC] = images_cropped[iIMC].astype(np.uint8)

    # Write to gif, simply add gif extension to folder to define gif path if the gif path is not specified
    if not output_gif_path:
        output_gif_path = folder_path + '.gif'
    imageio.mimsave(output_gif_path, images_cropped, fps=fps)

    # Optimize and overwrite gif
    optimize(output_gif_path)


# Function to compute volumes for image
def subjectPath2volumeTable(subject_path, print_bool=False):
    # Compute voxel numbers and volumes and output to table

    # Load image
    subject_image = nib.load(subject_path)
    subject = subject_image.get_fdata()

    # Get voxel volume
    voxel_volume = np.prod(subject_image.header['pixdim'][1:4])  # should be in mm^3

    # Calculate volumes
    [iVoxel, nVoxel] = np.unique(np.int64(np.round(subject)),
                                 return_counts=True)
    vVoxel = nVoxel * voxel_volume

    # # VTA check, remove later
    # VTA_volume = nVoxel[iVoxel == 1241]
    # print(f'subject path = {subject_path}, with VTA volume {VTA_volume}')

    # Output to DataFrame
    volume_table = pd.DataFrame(
        {'VolumeInteger': iVoxel,
         'VoxelNumber': nVoxel,
         'Volume': vVoxel})

    if print_bool:
        print(subject_image.header['pixdim'][1:4])
        print('voxel volume = ' + str(voxel_volume) + 'mm^3')
        print('total volume = ' + str(np.sum(vVoxel[iVoxel != 0])) + 'mm^3')

    return volume_table


# Function to compute volumes for image
def imageFLIRT2defField(image, flirt):
    M_image = image.affine[:3, :3]
    abc_image = image.affine[:3, 3]
    # dof = kwargs.get('dof', '12')
    # if dof == '6':
    #     print('rigid')
    #     M_flirt = flirt
    #     abc_flirt = 0
    # else:
    print('affine')
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
def zeroPadImage(input_3D_numpy, input_3D_numpy_template, padRatio):
    # crop
    nonzero_logical = np.where(input_3D_numpy_template > 0)
    edgeMax = np.max(nonzero_logical, axis=1)
    edgeMin = np.min(nonzero_logical, axis=1)
    edgeDis = edgeMax - edgeMin
    edgeDis_max = np.max(edgeDis)
    padLength = int(np.round(edgeDis_max * padRatio))

    print(f'edgeMax = {edgeMax}')
    print(f'edgeMin = {edgeMin}')
    print(f'edgeDis_max = {edgeDis_max}')
    print(f'padLength = {padLength}')

    output = input_3D_numpy[edgeMin[0]:edgeMax[0]+1, edgeMin[1]:edgeMax[1]+1, edgeMin[2]:edgeMax[2]+1]

    # check whether amount of nonzero elements is still equal after cropping
    print((f'Original number of nonzero elements = {np.sum(input_3D_numpy>0)}'))
    print((f'Template number of nonzero elements = {np.sum(input_3D_numpy_template>0)}'))
    print((f'Maximally cropped output number of nonzero elements = {np.sum(output>0)}'))

    print(f'Cropped input shape = {output.shape}')
    output = np.pad(output,
                    ((padLength, padLength),
                     (padLength, padLength),
                     (padLength, padLength)),
                    'constant')
    print((f'Zero padded output number of nonzero elements = {np.sum(output>0)}'))
    print(f'Zeropadded output shape = {output.shape}')

    # output_crop_index = edgeMin - padLength
    # print((f'Output crop index = {output_crop_index}'))
    #
    # input_crop_index = padLength - edgeMin
    # print((f'Input crop index = {input_crop_index}'))

    crop_index = (edgeMin, edgeMax, padLength)

    return output, crop_index

    # input_3D_numpy[crop_index[0]:, crop_index[1]:, crop_index[2]:] = output
    # output[crop_index[0]]

def unzeroPadImage(output, crop_index, input_3D_numpy):
    ## Convert possibly altered output back to input space by cropping and assigning to an input sized array

    # parse inputs
    edgeMin = crop_index[0]
    edgeMax = crop_index[1]
    padLength = crop_index[2]

    # remove zeropadding
    print(f'Output shape = {output.shape}')
    output = output[padLength:output.shape[0]-padLength,
                    padLength:output.shape[1]-padLength,
                    padLength:output.shape[2]-padLength]
    print(f'Unzeropadded output shape = {output.shape}')

    # Assign output to input sized array
    input_sized = np.zeros(input_3D_numpy.shape, dtype=input_3D_numpy.dtype)
    input_sized[edgeMin[0]:edgeMax[0]+1,
                edgeMin[1]:edgeMax[1]+1,
                edgeMin[2]:edgeMax[2]+1] = output
    print(f'Assigned input shape = {input_sized.shape}')

    return input_sized


def round_half_up(number, dec_places=0):
    sign = math.copysign(1, number)
    number_exact = abs(fractions.Fraction(number))
    shifted = number_exact * 10**dec_places
    shifted_trunc = int(shifted)
    if shifted - shifted_trunc >= fractions.Fraction(1, 2):
        result = (shifted_trunc + 1) / 10**dec_places
    else:
        result = shifted_trunc / 10**dec_places
    return sign * float(result)

def strPathList2List(strPathList):
    List = list(map(int, strPathList.strip('][').split(', ')))
    return List

def getChildStructures_mouse(group_ids_list, allen_table):
    # get all allen atlas children of group of ids
    # map(int, mouse_table_reference.loc[iVolume, 'structure_id_path_custom'].strip('][').split(', ')))[:-1]:
    child_list = list()
    for iRow in range(allen_table.shape[0]):
        structure_id_path_custom = list(map(int, allen_table.loc[iRow, 'structure_id_path_custom'].strip('][').split(', ')))
        structure_id_path_custom_in_group_ids = np.isin(structure_id_path_custom, group_ids_list)
        # print(structure_id_path_custom)
        if np.any(structure_id_path_custom_in_group_ids):
            first_child_id = np.where(structure_id_path_custom_in_group_ids)[0][0]
            # print(f'first_child_id = {first_child_id}')
            # print(f'structure_id_path_custom = {structure_id_path_custom}')
            append_this = np.array(structure_id_path_custom[first_child_id:])
            # print(f'append_this = {append_this}')
            append_this_len = len(append_this)
            # print(f'append_this_len = {append_this_len}')
            if append_this_len > 0:
                # print('APPENDING')
                child_list.append(append_this)
    child = np.unique(np.concatenate(child_list))
    return child

def map_annotation_fromto(annotation, map_from, map_to):
    annotation_remapped = np.round(annotation)  # always annotation so should never be non-integer
    annotation_remapped = annotation_remapped.astype(int) # always annotation so should never be non-integer

    # Remap map_from integers to map_to values
    annotation_remapped_shape = annotation_remapped.shape
    annotation_remapped = annotation_remapped.reshape(-1)
    annotation_remapped = npi.remap(annotation_remapped, map_from, map_to)
    annotation_remapped = annotation_remapped.reshape(annotation_remapped_shape)

    return(annotation_remapped)


def load_image(input_path, is_annotation=False):
    print(f'Loading {input_path}')
    input_image = nib.load(input_path)
    input = input_image.get_fdata()
    if is_annotation:
        input = np.round(input).astype(int)

@jit(nopython=False)
def computeJacobianDet(back_field, reference_annotation_in_group):
    jacobian_xx = np.gradient(np.squeeze(back_field[:, :, :, 0]), axis=0)
    jacobian_xy = np.gradient(np.squeeze(back_field[:, :, :, 0]), axis=1)
    jacobian_xz = np.gradient(np.squeeze(back_field[:, :, :, 0]), axis=2)
    jacobian_yx = np.gradient(np.squeeze(back_field[:, :, :, 1]), axis=0)
    jacobian_yy = np.gradient(np.squeeze(back_field[:, :, :, 1]), axis=1)
    jacobian_yz = np.gradient(np.squeeze(back_field[:, :, :, 1]), axis=2)
    jacobian_zx = np.gradient(np.squeeze(back_field[:, :, :, 2]), axis=0)
    jacobian_zy = np.gradient(np.squeeze(back_field[:, :, :, 2]), axis=1)
    jacobian_zz = np.gradient(np.squeeze(back_field[:, :, :, 2]), axis=2)
    jacobian_det = np.zeros(back_field.shape[0:3])
    for iX in range(back_field.shape[0]):
        # iX_percentage = (iX / back_field.shape[0]) * 100
        # print(f'iX_percentage = {iX_percentage}')
        for iY in range(back_field.shape[1]):
            for iZ in range(back_field.shape[2]):
                if reference_annotation_in_group[iX, iY, iZ]:
                    jacobian_det[iX, iY, iZ] = np.linalg.det(
                        np.array([[jacobian_xx[iX, iY, iZ], jacobian_xy[iX, iY, iZ], jacobian_xz[iX, iY, iZ]],
                                  [jacobian_yx[iX, iY, iZ], jacobian_yy[iX, iY, iZ], jacobian_yz[iX, iY, iZ]],
                                  [jacobian_zx[iX, iY, iZ], jacobian_zy[iX, iY, iZ], jacobian_zz[iX, iY, iZ]]]))
    return jacobian_det

def computeJacobianDet_nojit(back_field, reference_annotation_in_group):
    jacobian_xx = np.gradient(np.squeeze(back_field[:, :, :, 0]), axis=0)
    jacobian_xy = np.gradient(np.squeeze(back_field[:, :, :, 0]), axis=1)
    jacobian_xz = np.gradient(np.squeeze(back_field[:, :, :, 0]), axis=2)
    jacobian_yx = np.gradient(np.squeeze(back_field[:, :, :, 1]), axis=0)
    jacobian_yy = np.gradient(np.squeeze(back_field[:, :, :, 1]), axis=1)
    jacobian_yz = np.gradient(np.squeeze(back_field[:, :, :, 1]), axis=2)
    jacobian_zx = np.gradient(np.squeeze(back_field[:, :, :, 2]), axis=0)
    jacobian_zy = np.gradient(np.squeeze(back_field[:, :, :, 2]), axis=1)
    jacobian_zz = np.gradient(np.squeeze(back_field[:, :, :, 2]), axis=2)
    jacobian_det = np.zeros(back_field.shape[0:3])
    for iX in range(back_field.shape[0]):
        # iX_percentage = (iX / back_field.shape[0]) * 100
        # print(f'iX_percentage = {iX_percentage}')
        for iY in range(back_field.shape[1]):
            for iZ in range(back_field.shape[2]):
                if reference_annotation_in_group[iX, iY, iZ]:
                    jacobian_det[iX, iY, iZ] = np.linalg.det(
                        np.array([[jacobian_xx[iX, iY, iZ], jacobian_xy[iX, iY, iZ], jacobian_xz[iX, iY, iZ]],
                                  [jacobian_yx[iX, iY, iZ], jacobian_yy[iX, iY, iZ], jacobian_yz[iX, iY, iZ]],
                                  [jacobian_zx[iX, iY, iZ], jacobian_zy[iX, iY, iZ], jacobian_zz[iX, iY, iZ]]]))
    return jacobian_det

def computeCohenD_WT_KO(nWT, nKO, std_WT, std_KO, mean_WT, mean_KO, nIterBootstrap=0):
    """
    Calculate CohenD effect size with pooled variance and hedge correction.
    """
    N = nWT + nKO
    hedge_correction = (N - 3) / (N - 2.25)
    S_p = np.sqrt(((nWT - 1) * np.power(std_WT, 2) + (nKO - 1) * np.power(std_KO, 2)) / (N - 2))

    cohenD = ((mean_WT - mean_KO) / S_p) * hedge_correction

    if nIterBootstrap != 0:
        cohenD_BS = np.empty(nIterBootstrap)
        for iBS in range(nIterBootstrap):
            WT_BS = np.random.normal(mean_WT, std_WT, nWT)
            KO_BS = np.random.normal(mean_KO, std_KO, nKO)

            mean_WT_BS = np.mean(WT_BS)
            mean_KO_BS = np.mean(KO_BS)

            std_WT_BS = np.std(mean_WT_BS)
            std_KO_BS = np.std(mean_KO_BS)

            S_p = np.sqrt(((nWT - 1) * np.power(std_WT_BS, 2) + (nKO - 1) * np.power(std_KO_BS, 2)) / (nKO + nWT - 2))
            cohenD_BS[iBS] = ((mean_WT_BS - mean_KO_BS) / S_p) * hedge_correction # positve is greater WT mean, negative is smaller WT

        cohenD_CI = [np.quantile(cohenD_BS, .025), np.quantile(cohenD_BS, .975)]

        return cohenD, cohenD_CI

    else:

        return cohenD

def voxCoords2fslCoords(voxCoords, image, dims):
    # If the voxel to mm mapping has a positive determinant convert x dim according to fsl specifications
    if np.linalg.det(image.get_qform()) > 0:
        voxCoords[0] = image.shape[0] - 1 - voxCoords[0]

    fslCoords = voxCoords
    fslCoords[0:3] = voxCoords[0:3] * dims

    return fslCoords

def fslCoords2voxCoords(fslCoords, image, dims):

    voxCoords = fslCoords
    voxCoords[0:3] = fslCoords[0:3] / dims

    # If the voxel to mm mapping has a positive determinant convert x dim according to fsl specifications
    if np.linalg.det(image.get_qform()) > 0:
        voxCoords[0] = image.shape[0] - 1 - voxCoords[0]

    return voxCoords

def timeFunc(Func, input_list, nIter):
    elapsed_time = list()
    for iIter in range(nIter):
        t = time.process_time()

        Func(*input_list)

        elapsed_time.append(time.process_time() - t)

    print(f'elapsed_time = {np.mean(elapsed_time)}')
    return elapsed_time

def CrawfordHowell(case, control):
    tval = (case - np.mean(control)) / (np.std(control) * np.sqrt((len(control) + 1) / len(control)))
    degfree = len(control) - 1
    pval = 2 * stats.t.sf(np.abs(tval), df=degfree)  # two-tailed p-value
    return tval, degfree, pval

