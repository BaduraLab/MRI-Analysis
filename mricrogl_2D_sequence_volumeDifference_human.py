import gl
import sys
import os
print(sys.version)
print(gl.version())
gl.resetdefaults()

help(gl)

# Load main image and overlay
mep_path = os.path.join('/', 'mnt', 'tosh', 'Projects', 'MEP', 'mep-scripts')
reference_path = os.path.join(mep_path, 'Data', 'Human', 'Reference')
analysis_path = os.path.join(mep_path, 'Data', 'Human', 'Analysis')
main_path = os.path.join(reference_path, 'standard', 'MNI152_T1_0.5mm.nii.gz')
overlay1_path = os.path.join(analysis_path, 'Lobules-SUIT_aFD_volIncrease.nii.gz')
overlay2_path = os.path.join(analysis_path, 'Lobules-SUIT_aFD_volDecrease.nii.gz')
overlay3_path = os.path.join(analysis_path, 'prob_atlas_bilateral_thrarg_0_aFD_volIncrease.nii.gz')
overlay4_path = os.path.join(analysis_path, 'prob_atlas_bilateral_thrarg_0_aFD_volDecrease.nii.gz')
gl.loadimage(main_path)
gl.overlayload(overlay1_path)
gl.overlayload(overlay2_path)
gl.overlayload(overlay3_path)
gl.overlayload(overlay4_path)
gl.opacity(1, 50)
gl.opacity(2, 50)
gl.opacity(3, 50)
gl.opacity(4, 50)
gl.minmax(1, 0, 0.25)
gl.minmax(2, 0, 0.25)
gl.minmax(3, 0, 0.25)
gl.minmax(4, 0, 0.25)
gl.colorname(1, '1red')
gl.colorname(2, '3blue')
gl.colorname(3, '1red')
gl.colorname(4, '3blue')
output_dir_path = os.path.join(analysis_path, 'imageSequenceFolders', 'aFD_InDe')
if not os.path.exists(output_dir_path):
	os.mkdir(output_dir_path)

y_min = -120
y_max = 90
y_step = 1

y = y_min
count = 0
while y <= y_max:
	print(y)
	count = count + 1
	print(count)

	gl.orthoviewmm(-3, y, -14)
	gl.view(2)
	gl.linewidth(0)
	gl.colorbarposition(0)

	filepath = os.path.join(output_dir_path, 'imagenew' + '_' + str(round(count)).rjust(3, '0')) + '_' + str(y)
	gl.savebmp(filepath)

	y = y + y_step
	


