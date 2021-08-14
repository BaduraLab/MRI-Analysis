# For a given image, overlays and reference voxel coordinates, take all relevant images (mostly masks and background)
# It might be a time saver and give improved image quality to split images on axial-sagittal-coronal
import gl
import sys
import os
import glob
print(sys.version)
print(gl.version())
# gl.resetdefaults()
import pandas as pd
import ast

sysBase_path = os.path.join('/', 'mnt', 'tosh', 'Projects', 'MEP', 'mep-scripts')
data_path = os.path.join(sysBase_path, 'Data')
analysis_path = os.path.join(data_path, 'Analysis')
coords_table_path = os.path.join(analysis_path, 'coords_table.csv')
coords_table = pd.read_csv(coords_table_path, converters={"ref_Coords": ast.literal_eval,
                                                          "flirtedRigid_Coords": ast.literal_eval,
                                                          "native_Coords": ast.literal_eval})
data_human_path = os.path.join(data_path, 'Human', 'Processed')
data_mouse_path = os.path.join(data_path, 'Mouse', 'Processed_Old')
# data_path = os.path.join('/', 'mnt', 'tosh', 'Projects', 'MEP', 'mep-scripts', 'Data', 'Mouse', 'Processed_Old')

# Coords table filtering
coords_table = coords_table[[strID in ['human_patient_SN_CE', 'human_control_SN_CE'] for strID in coords_table['strID']]]

# Load table and loop through rows, no data_subject_path or subject_folder list necessary
for iRow in range(coords_table.shape[0]):
    coords_table_row = coords_table.iloc[iRow]
    data_subject_path = coords_table_row.loc['data_subject_path']
    subject = data_subject_path.split(os.sep)[-1]
    species = data_subject_path.split(os.sep)[1]

    output_folder = os.path.join(analysis_path, coords_table_row.loc['strID'])
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
# for data_subject_path in [data_human_path, data_mouse_path]:
#     subject_folder_list = glob.glob(os.path.join(data_subject_path, '*'))


    transform_level_list = ['native', 'native-skull', 'flirtedRigid', 'flirted', 'synned', 'native-adjusted']
    for transform_level in transform_level_list:
        print(f'transform_level = {transform_level}')
        # for subject_folder in subject_folder_list:
        #     subject = subject_folder.split(os.sep)[-1]

        output_transform_folder = os.path.join(output_folder, transform_level)
        if not os.path.exists(output_transform_folder):
            os.makedirs(output_transform_folder)

        # with or without skull? skull removal images of course with skull

        # Get template which only differ per transform level
        if transform_level == 'native':
            if species == 'Human':
                template_path = os.path.join(sysBase_path, data_subject_path, subject + '_reoriented.nii.gz')
            else:
                template_path = os.path.join(sysBase_path, data_subject_path, 'FLIRT', subject + '_inmasked.nii.gz')
            Coords = coords_table_row.loc['native_Coords']
            # Coords = eval(coords_table_row.loc['native_Coords'])
            # print(f'Coords_str = {Coords_str}')
            print(f'Coords = {Coords}')
        if transform_level == 'native-skull':
            if species == 'Human':
                template_path = os.path.join(sysBase_path, data_subject_path, subject + '_skull_reoriented.nii.gz')
            else:
                template_path = os.path.join(sysBase_path, data_subject_path, 'FLIRT', subject + '.nii.gz')
            Coords = coords_table_row.loc['native_Coords']
        elif transform_level == 'flirtedRigid':
            template_path = os.path.join(sysBase_path, data_subject_path, 'retransform', subject + '_' + transform_level + '.nii.gz')
            Coords = coords_table_row.loc['flirtedRigid_Coords']
        elif transform_level == 'flirted':
            template_path = os.path.join(sysBase_path, data_subject_path, 'retransform', subject + '_' + transform_level + '.nii.gz')
            Coords = coords_table_row.loc['ref_Coords']
        elif transform_level == 'synned':
            template_path = os.path.join(sysBase_path, data_subject_path, 'retransform', subject + '_' + transform_level + '.nii.gz')
            Coords = coords_table_row.loc['ref_Coords']
        else: # native-adjusted
            if species == 'Human':
                template_path = os.path.join(sysBase_path, data_subject_path, subject + '_reoriented.nii.gz')
            else:
                template_path = os.path.join(sysBase_path, data_subject_path, 'FLIRT', subject + '_inmasked.nii.gz')
            Coords = coords_table_row.loc['native_Coords']

        # Get annotations
        data_Isolations_path = os.path.join(sysBase_path, data_subject_path, 'Isolations_' + transform_level)
        data_Isolations_annotation_path_list = glob.glob(os.path.join(data_Isolations_path, '*', '*annotation*.nii.gz'))

        # Define images to load and unload
        data_load_path_list = [template_path] + data_Isolations_annotation_path_list

        for data_load_path in data_load_path_list:

            gl.overlayloadsmooth(0)
            gl.loadimage(data_load_path)
            gl.overlayloadsmooth(0)

            # overlay_path_list = glob.glob(os.path.join(data_structureGroup_path,
            #                        subject + '_annotation_*' + structureAcronym + '.nii.gz'))
            # gl.overlayload(overlay_path_list[0])
            # gl.opacity(1, 100)

            if data_load_path == template_path:
                structureNameID = 'background'
                help(gl.minmax)
                if species == 'Mouse':
                    gl.minmax(0, 0, 130)
                    # gl.minmax(0, 0, 0.2)
            else:
                structureNameID = data_load_path.split(os.sep)[-1].split('_')[-1]
                structureNameID_addition = data_load_path.split(os.sep)[-1].split('_')[-2]
                if (structureNameID_addition != transform_level) & (structureNameID_addition != 'thrarg'):
                    structureNameID = structureNameID_addition + '_' + structureNameID
                    structureNameID_addition = data_load_path.split(os.sep)[-1].split('_')[-3]
                    if (structureNameID_addition != transform_level) & (structureNameID_addition != 'thrarg'):
                        structureNameID = structureNameID_addition + '_' + structureNameID
                gl.minmax(0, 0, 0.1)

            # gl.colorname(1, 'FABulous_PAT')

            # loop over viewVals
            for viewVal in [1, 2, 4]:
                print(f'viewVal={viewVal}')
                if viewVal == 1:
                    print(f'Coords[0] = {Coords[0]}')
                    print(f'Coords[1] = {Coords[1]}')
                    print(f'Coords[2] = {Coords[2]}')
                    print(f'z_center = [{Coords[0]}, {Coords[1]}, {Coords[2]}]')
                    gl.orthoviewmm(Coords[0], Coords[1], Coords[2])
                    gl.view(viewVal)
                    gl.linewidth(0)
                    gl.colorbarposition(0)
                    # gl.minmax(0, z_minmax[0], z_minmax[1])
                    viewString = 'axial'
                elif viewVal == 2:
                    print(f'Coords[0] = {Coords[0]}')
                    print(f'Coords[1] = {Coords[1]}')
                    print(f'Coords[2] = {Coords[2]}')
                    print(f'y_center = [{Coords[0]}, {Coords[1]}, {Coords[2]}]')
                    gl.orthoviewmm(Coords[0], Coords[1], Coords[2])
                    gl.view(viewVal)
                    gl.linewidth(0)
                    gl.colorbarposition(0)
                    # gl.minmax(0, y_minmax[0], y_minmax[1])
                    viewString = 'coronal'
                elif viewVal == 4:
                    print(f'Coords[0] = {Coords[0]}')
                    print(f'Coords[1] = {Coords[1]}')
                    print(f'Coords[2] = {Coords[2]}')
                    print(f'x_center = [{Coords[0]}, {Coords[1]}, {Coords[2]}]')
                    gl.orthoviewmm(Coords[0], Coords[1], Coords[2])
                    gl.view(viewVal)
                    gl.linewidth(0)
                    gl.colorbarposition(0)
                    # gl.minmax(0, x_minmax[0], x_minmax[1])
                    viewString = 'sagittal'

                # gl.opacity(1, 100)
                # gl.bmpzoom(2)
                filepath = os.path.join(output_transform_folder, viewString + '_' + structureNameID)
                print(f'filepath = {filepath}')
                gl.savebmp(filepath)


        # cross section full screen possibly zoomed in with overlay with around 0.4 opacity

        # rigid human
        # structure isolation human edit
        # cross human script

        # 3D human/mouse already done I think vol script with rotations
        # just choose orientations for custom pictures, don't worry on automating because view is more important

