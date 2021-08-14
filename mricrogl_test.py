import gl
import sys
print(sys.version)
print(gl.version())
#gl.resetdefaults()
import os
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


gl.overlayloadsmooth(0)
gl.loadimage(os.path.join(data_mouse_path, 'KO_3B_9', 'KO_3B_9.nii.gz'))
gl.overlayload(os.path.join(data_mouse_path, 'KO_3B_9', 'allen_annotation_invsynned_to_KO_3B_9.nii.gz'))


gl.overlayloadsmooth(0)
