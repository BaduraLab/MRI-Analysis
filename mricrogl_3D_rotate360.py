# Loop through rotations and save frame
import gl
gl.view(64)

import sys
from datetime import datetime
import os
strtime = str(datetime.now())
nRot = 100
output_path = os.path.join('/', 'mnt', 'tosh', 'Projects', 'MEP', 'mep-scripts', 'Data', 'Analysis',
                           'imageSequenceFolders', 'rotate_3D_' + strtime)
if not os.path.exists(output_path):
    os.mkdir(output_path)

gl.shadername('default')

for iRot in range(nRot):
    print(iRot)
    Rot = iRot * (360 / nRot)
    print(Rot)
    gl.azimuthelevation(round(Rot), 20)
    filepath = output_path + os.sep + 'image_' + str(round(Rot)).rjust(3, '0')
    gl.savebmp(filepath)
    print(filepath)
