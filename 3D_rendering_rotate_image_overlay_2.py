import gl
import sys
print(sys.version)
print(gl.version())
import os

help(gl)
gl.linewidth(1)



gl.view(64)

round(0.1)
nRot = 50
data_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Data/Mouse/Processed_New'
subject = 'KO_6'
input_path = data_path + '/'+subject+'/original_nonisolated_cerebellum.nii.gz'
overlay_path = data_path + '/'+subject+'/annotation_nonisolated_cerebellum.nii.gz'
#overlay_path = data_path + '/'+subject+'/allen_annotation_invsynned_to_'+subject+'_adjusted_cerebellum.nii.gz'
output_path = 'C:/Users/enzo/Downloads/Study/Current Courses/MEP/Pictures/'+subject+'_50rot_overlay_cutout'
if not os.path.exists(output_path):
	os.mkdir(output_path)
gl.loadimage(input_path)
gl.overlayload(overlay_path)
gl.opacity(1, 40)
gl.colorname(1, 'x_rain')
#gl.fullscreen(0)
gl.cutout(0.5, 0.5, 0.5, 0, 0, 1)

help(gl.colorname)

for iRot in range(nRot):
	print(iRot)
	Rot = iRot * (360/nRot)
	print(Rot)
	gl.azimuthelevation(round(Rot), 20)
	gl.savebmp(output_path+'/imagenew'+str(round(Rot)))






# gl.resetdefaults()
# gl.viewcoronal(True)
#gl.view(2)
#gl.orthoviewmm(0.5,0.5,0.5)