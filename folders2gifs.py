from functions import imageFolder2gif
import glob
import os
import PIL

# save all image sequences saved in analysis image sequence folder
data_path = os.path.join('Data', 'Mouse', 'Processed_Old')
analysis_path = os.path.join('Data', 'Mouse', 'Analysis', 'imageSequenceFolders')

gif_duration = 300

for folderpath in list(set(glob.glob(os.path.join(analysis_path, '*')))-set(glob.glob(os.path.join(analysis_path, '*.gif')))):
    print(folderpath)

    filepath_list = glob.glob(os.path.join(folderpath, '*'))
    im_list = [PIL.Image.open(filepath) for filepath in filepath_list]
    im_list[0].save(folderpath+'.gif',
                    save_all=True,
                    append_images=im_list[1:],
                    optimize=False, duration=gif_duration, loop=0)

# Do the same for human data
data_path = os.path.join('Data', 'Human', 'Processed')
analysis_path = os.path.join('Data', 'Human', 'Analysis', 'imageSequenceFolders')

for folderpath in list(set(glob.glob(os.path.join(analysis_path, '*')))-set(glob.glob(os.path.join(analysis_path, '*.gif')))):
    print(folderpath)

    filepath_list = glob.glob(os.path.join(folderpath, '*'))
    im_list = [PIL.Image.open(filepath) for filepath in filepath_list]
    im_list[0].save(folderpath + '.gif',
                    save_all=True,
                    append_images=im_list[1:],
                    optimize=False, duration=gif_duration, loop=0)
