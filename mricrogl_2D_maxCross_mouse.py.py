import gl
import sys
import os
import glob
print(sys.version)
print(gl.version())
gl.resetdefaults()
import pandas as pd
import ast

data_path = os.path.join('/', 'mnt', 'tosh', 'Projects', 'MEP', 'mep-scripts', 'Data', 'Mouse', 'Processed_Old')
subject_folder_list = glob.glob(os.path.join(data_path, '*'))

group_name_list = ['ci', 'expSel']

for group_name in group_name_list:
    for subject_folder in subject_folder_list:
        subject = subject_folder.split(os.sep)[-1]

        data_structureGroup_path = os.path.join(subject_folder, 'Isolations', group_name)
        expSel_table_path = os.path.join(data_structureGroup_path, group_name + '.csv')
        expSel_table = pd.read_csv(expSel_table_path)

        for iRow in range(expSel_table.shape[0]):
            structureAcronym = expSel_table.iloc[iRow]['acronym']

            gl.loadimage(os.path.join(subject_folder, 'FLIRT', subject + '_inmasked.nii.gz'))
            gl.overlayloadsmooth(0)
            gl.overlayload(os.path.join(data_structureGroup_path, 'allen_annotation_invsynned_to_' + subject + '_adjusted_cerebellum_lobular_' + structureAcronym + '.nii.gz'))
            # gl.opacity(1, 100)
            gl.minmax(1, 0, 0.2)
            if subject[0:2] == 'WT':
                gl.colorname(1, 'FABulous_WT')
            else:
                gl.colorname(1, 'FABulous_PAT')

            x_center = ast.literal_eval(expSel_table.iloc[iRow]['area_x_center_Coord'])
            y_center = ast.literal_eval(expSel_table.iloc[iRow]['area_y_center_Coord'])
            z_center = ast.literal_eval(expSel_table.iloc[iRow]['area_z_center_Coord'])

            x_minmax = ast.literal_eval(expSel_table.iloc[iRow]['template_x_minmax'])
            y_minmax = ast.literal_eval(expSel_table.iloc[iRow]['template_y_minmax'])
            z_minmax = ast.literal_eval(expSel_table.iloc[iRow]['template_z_minmax'])

            # loop over viewVals
            gl.opacity(1, 100)
            for viewVal in [1, 2, 4]:
                print(f'viewVal={viewVal}')
                if viewVal == 1:
                    gl.orthoviewmm(z_center[0], z_center[1], z_center[2])
                    gl.view(viewVal)
                    gl.linewidth(0)
                    gl.colorbarposition(0)
                    gl.minmax(0, z_minmax[0], z_minmax[1])
                    viewString = 'axial'
                elif viewVal == 2:
                    gl.orthoviewmm(y_center[0], y_center[1], y_center[2])
                    gl.view(viewVal)
                    gl.linewidth(0)
                    gl.colorbarposition(0)
                    gl.minmax(0, y_minmax[0], y_minmax[1])
                    viewString = 'coronal'
                elif viewVal == 4:
                    gl.orthoviewmm(x_center[0], x_center[1], x_center[2])
                    gl.view(viewVal)
                    gl.linewidth(0)
                    gl.colorbarposition(0)
                    gl.minmax(0, x_minmax[0], x_minmax[1])
                    viewString = 'sagittal'

                # gl.bmpzoom(2)
                filepath = os.path.join(data_structureGroup_path, 'image_' + viewString + '_' + structureAcronym)
                print(f'filepath = {filepath}')
                gl.savebmp(filepath)
                # gl.orthoviewmm()

            # loop over viewVals
            gl.opacity(1, 0)
            for viewVal in [1, 2, 4]:
                print(f'viewVal={viewVal}')
                if viewVal == 1:
                    gl.orthoviewmm(z_center[0], z_center[1], z_center[2])
                    gl.view(viewVal)
                    gl.linewidth(0)
                    gl.colorbarposition(0)
                    gl.minmax(0, z_minmax[0], z_minmax[1])
                    viewString = 'axial'
                elif viewVal == 2:
                    gl.orthoviewmm(y_center[0], y_center[1], y_center[2])
                    gl.view(viewVal)
                    gl.linewidth(0)
                    gl.colorbarposition(0)
                    gl.minmax(0, y_minmax[0], y_minmax[1])
                    viewString = 'coronal'
                elif viewVal == 4:
                    gl.orthoviewmm(x_center[0], x_center[1], x_center[2])
                    gl.view(viewVal)
                    gl.linewidth(0)
                    gl.colorbarposition(0)
                    gl.minmax(0, x_minmax[0], x_minmax[1])
                    viewString = 'sagittal'

                # gl.bmpzoom(2)
                filepath = os.path.join(data_structureGroup_path, 'image_' + viewString + '_' + structureAcronym + '_transparent')
                print(f'filepath = {filepath}')
                gl.savebmp(filepath)


        # cross section full screen possibly zoomed in with overlay with around 0.4 opacity

        # rigid human
        # structure isolation human edit
        # cross human script

        # 3D human/mouse already done I think vol script with rotations
        # just choose orientations for custom pictures, don't worry on automating because view is more important

