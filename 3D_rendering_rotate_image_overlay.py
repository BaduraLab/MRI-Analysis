import gl
import sys
print(sys.version)
print(gl.version())
import os

help(gl)
gl.linewidth(1)



gl.view(64)

round(0.1)
nRot = 25
data_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse/Processed_New'
input_path = data_path + '/WT_50/FLIRT/WT_50_inmasked.nii.gz'
overlay_path = data_path + '/WT_50/allen_annotation_invsynned_to_WT_50_adjusted_cerebellum.nii.gz'
output_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Pictures/Rotation'
if not os.path.exists(output_path):
	os.mkdir(output_path)
gl.loadimage(input_path)
gl.overlayload(overlay_path)
gl.opacity(1, 40)
gl.colorname(1, 'x_rain')
#gl.fullscreen(0)

for iRot in range(nRot):
	print(iRot)
	Rot = iRot * (360/nRot)
	print(Rot)
	gl.azimuthelevation(round(Rot), 20)
	gl.savebmp(output_path+'/doesthisworknow'+str(round(Rot)))






# gl.resetdefaults()
# gl.viewcoronal(True)
#gl.view(2)
#gl.orthoviewmm(0.5,0.5,0.5)