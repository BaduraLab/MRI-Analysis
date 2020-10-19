import gl
import sys
import os
import glob
from pathlib import Path
print(sys.version)
print(gl.version())
gl.resetdefaults()

# help(gl)

# Load main image and overlay
mep_path = os.path.join('D:'+os.sep, 'Users', 'enzo', 'Documents', 'Projects', 'MEP', 'mep-scripts')
reference_path = os.path.join(mep_path, 'Data', 'Human', 'Reference')
analysis_path = os.path.join(mep_path, 'Data', 'Human', 'Analysis')
data_path = os.path.join(mep_path, 'Data', 'Human', 'Processed')
subject_path_list = glob.glob(os.path.join(data_path, '*'))

print(subject_path_list)

for subject_path in subject_path_list:
	subject = subject_path.split(os.sep)[-1]

	main_path = os.path.join(subject_path, subject+'_reoriented.nii.gz')
	overlay1_path = os.path.join(subject_path, subject+'_annotation_orsuit_thrarg_adjusted_lobular_mc.nii.gz')

	gl.loadimage(main_path)
	gl.overlayloadsmooth(0)
	gl.overlayload(overlay1_path)

	gl.opacity(1, 50)
	gl.minmax(1, 0, 34)
	gl.colorname(1, 'x_rain')

	output_dir_path = os.path.join(analysis_path, 'imageSequenceFolders', subject+'_orthoview')
	if not os.path.exists(output_dir_path):
		os.mkdir(output_dir_path)

	if subject == 'patient':
		y_min = -84
		y_max = -14
	else:
		y_min = -80
		y_max = -10
	
	y_step = 1

	y = y_min
	count = 0
	while y <= y_max:
		print(y)
		count = count + 1
		print(count)

		gl.orthoviewmm(0, y, 0) # why at first (-3, y, -14)?
		gl.view(2)
		gl.linewidth(0)
		gl.colorbarposition(0)

		filepath = os.path.join(output_dir_path, 'imagenew' + '_' + str(round(count)).rjust(3, '0')) + '_' + str(y)
		gl.savebmp(filepath)

		y = y + y_step
