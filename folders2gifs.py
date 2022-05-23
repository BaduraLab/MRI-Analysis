#!/usr/bin/python3
"""  Convert folders of images to GIFs
"""

__author__ = "Enzo Nio"
__version__ = "1.0.0"
__maintainer__ = "Enzo Nio"

from functions import imageFolder2gif
import glob
import os
import numpy as np
import imageio

# logical whether to run gif creation for sets of data
run_logical = [False, False, False, True, False] # mouse, human and general analyses respectively
gif_duration = 300
analysis_path_list = [os.path.join('Data', 'Mouse', 'Analysis', 'imageSequenceFolders'),
                      os.path.join('Data', 'Human', 'Analysis', 'imageSequenceFolders'),
                      os.path.join('Data', 'Human', 'Analysis', 'imageSequenceFolders_isolations'),
                      os.path.join('Data', 'Analysis', 'imageSequenceFolders', 'runFolders'),
                      os.path.join('Data', 'Analysis', 'imageSequenceFolders')]

# Go through analysis folders and convert image sequence folders into gifs
for iAnalysis in range(len(analysis_path_list)):
    if run_logical[iAnalysis]:
        # Get analysis path
        analysis_path = analysis_path_list[iAnalysis]

        # List image sequence folders, ignore previously created gif of image sequence folders and overwrite these gifs
        for folderpath in list(set(glob.glob(os.path.join(analysis_path, '*')))-set(glob.glob(os.path.join(analysis_path, '*.gif')))):
            print(folderpath)
            # Convert image sequence folder to gif
            imageFolder2gif(folderpath)



# if run_logical[1]:
#     # Do the same for human data
#     analysis_path = os.path.join('Data', 'Human', 'Analysis', 'imageSequenceFolders')
#
#     # List image sequence folders, ignore previously created gif of image sequence folders and overwrite these gifs
#     for folderpath in list(set(glob.glob(os.path.join(analysis_path, '*')))-set(glob.glob(os.path.join(analysis_path, '*.gif')))):
#         print(folderpath)
#
#         # List image paths in image sequence folder
#         filepath_list = glob.glob(os.path.join(folderpath, '*'))
#
#         # Read sorted images
#         images = []
#         ymin = []; xmin = []
#         ymax = []; xmax = []
#         sorted_degrees_indices = np.argsort([int(filepath.split('.')[-2].split('_')[-1]) for filepath in filepath_list])
#         for iFilename in sorted_degrees_indices:
#             filename = filepath_list[iFilename]
#             print(filename)
#             images.append(imageio.imread(filename))
#             nonzero_yxCoord = np.where(np.max(images[-1], 2))
#
#             # Add minimum and maximum y and x coordinates for cropping
#             ymin.append(np.min(nonzero_yxCoord[0]))
#             ymax.append(np.max(nonzero_yxCoord[0]))
#             xmin.append(np.min(nonzero_yxCoord[1]))
#             xmax.append(np.max(nonzero_yxCoord[1]))
#
#         # Calculate outermost cropping box corners
#         ymin = np.min(ymin)
#         ymax = np.max(ymax) + 1
#         xmin = np.min(xmin)
#         xmax = np.max(xmax) + 1
#
#         # Crop images
#         images_cropped = []
#         for im in images:
#             images_cropped.append(im[ymin:ymax, xmin:xmax])
#
#         # Write to gif
#         gif_path = folderpath + '.gif'
#         imageio.mimsave(gif_path, images_cropped)
#
#         # Optimize and overwrite gif
#         optimize(gif_path)
#
#
#
# if run_logical[2]:
#     # Do the same for general data
#     analysis_path = os.path.join('Data', 'Analysis', 'imageSequenceFolders')
#
#     # List image sequence folders, ignore previously created gif of image sequence folders and overwrite these gifs
#     for folderpath in list(set(glob.glob(os.path.join(analysis_path, '*')))-set(glob.glob(os.path.join(analysis_path, '*.gif')))):
#         print(folderpath)
#
#         # List image paths in image sequence folder
#         filepath_list = glob.glob(os.path.join(folderpath, '*'))
#
#         # Read sorted images
#         images = []
#         ymin = []; xmin = []
#         ymax = []; xmax = []
#         sorted_degrees_indices = np.argsort([int(filepath.split('.')[-2].split('_')[-1]) for filepath in filepath_list])
#         for iFilename in sorted_degrees_indices:
#             filename = filepath_list[iFilename]
#             print(filename)
#             images.append(imageio.imread(filename))
#             nonzero_yxCoord = np.where(np.max(images[-1], 2))
#
#             # Add minimum and maximum y and x coordinates for cropping
#             ymin.append(np.min(nonzero_yxCoord[0]))
#             ymax.append(np.max(nonzero_yxCoord[0]))
#             xmin.append(np.min(nonzero_yxCoord[1]))
#             xmax.append(np.max(nonzero_yxCoord[1]))
#
#         # Calculate outermost cropping box corners
#         ymin = np.min(ymin)
#         ymax = np.max(ymax) + 1
#         xmin = np.min(xmin)
#         xmax = np.max(xmax) + 1
#
#         # Crop images
#         images_cropped = []
#         for im in images:
#             images_cropped.append(im[ymin:ymax, xmin:xmax])
#
#         # Write to gif
#         gif_path = folderpath + '.gif'
#         imageio.mimsave(gif_path, images_cropped)
#
#         # Optimize and overwrite gif
#         optimize(gif_path)



# data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
# # Do the same for human data
# data_path = os.path.join('Data', 'Human', 'Processed')
#
# for folderpath in list(set(glob.glob(os.path.join(analysis_path, '*')))-set(glob.glob(os.path.join(analysis_path, '*.gif')))):
#     print(folderpath)
#
#     filepath_list = glob.glob(os.path.join(folderpath, '*'))
#     im_list = [PIL.Image.open(filepath) for filepath in filepath_list]
#     im_list[0].save(folderpath + '.gif',
#                     save_all=True,
#                     append_images=im_list[1:],
#                     optimize=True, duration=gif_duration, loop=0)

# sequence = []
#
# im = Image.open(....)
#
# # im is your original image
# frames = [frame.copy() for frame in ImageSequence.Iterator(im)]
#
# # write GIF animation
# fp = open("out.gif", "wb")
# gifmaker.makedelta(fp, frames)
# fp.close()

# from PIL import ImageSequence
# from PIL import Image
# from pillow import gifmaker
# import PIL
