import gl
import sys
import os
print(sys.version)
print(gl.version())
gl.resetdefaults()

help(gl)

# Load main image and overlay
mep_path = os.path.join('C:/', 'Users', 'Enzon', 'Documents', 'Projects', 'MEP', 'mep-scripts')
reference_path = os.path.join(mep_path, 'Data', 'Human', 'Reference')
analysis_path = os.path.join(mep_path, 'Data', 'Human', 'Analysis')
main_path = os.path.join(reference_path, 'average_template_50_reoriented.nii.gz')
overlay1_path = os.path.join(analysis_path, 'annotation_50_reoriented_pVal_inv_sig_volIncrease.nii.gz')
overlay2_path = os.path.join(analysis_path, 'annotation_50_reoriented_pVal_inv_sig_volDecrease.nii.gz')
gl.loadimage(main_path)
gl.overlayload(overlay1_path)
gl.overlayload(overlay2_path)
gl.opacity(1, 50)
gl.opacity(2, 50)
gl.minmax(1, 0, 3)
gl.minmax(2, 0, 3)
gl.colorname(1, '3blue')
gl.colorname(2, '1red')
output_dir_path = os.path.join(analysis_path, 'imageSequenceFolders', 'pVal_inv_sig_InDe')
if not os.path.exists(output_dir_path):
	os.mkdir(output_dir_path)

y_min = -13.15
y_max = 0
y_step = 0.05

y = y_min
while y <= y_max:
	print(y)

	gl.orthoviewmm(0, y, -4)
	gl.view(2)
	gl.linewidth(0)

	filepath = os.path.join(output_dir_path, 'imagenew'+'_'+str(round(y*100)).rjust(4, '0'))
	gl.savebmp(filepath)

	y = y + y_step
	


