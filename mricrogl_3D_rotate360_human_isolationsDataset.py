import gl
import sys
print(sys.version)
print(gl.version())
import os
import glob

#sys.path.append(['C:\\Users\\Enzon\\Documents\\Projects\\MEP\\mep-scripts', 'C:\\Program Files\\JetBrains\\PyCharm 2020.1.2\\plugins\\python\\helpers\\pydev', 'C:\\Users\\Enzon\\Documents\\Projects\\MEP\\mep-scripts', 'C:\\Program Files\\JetBrains\\PyCharm 2020.1.2\\plugins\\python\\helpers\\pycharm_display', 'C:\\Program Files\\JetBrains\\PyCharm 2020.1.2\\plugins\\python\\helpers\\third_party\\thriftpy', 'C:\\Program Files\\JetBrains\\PyCharm 2020.1.2\\plugins\\python\\helpers\\pydev', 'C:\\Users\\Enzon\\AppData\\Local\\Programs\\Python\\Python37\\python37.zip', 'C:\\Users\\Enzon\\AppData\\Local\\Programs\\Python\\Python37\\DLLs', 'C:\\Users\\Enzon\\AppData\\Local\\Programs\\Python\\Python37\\lib', 'C:\\Users\\Enzon\\AppData\\Local\\Programs\\Python\\Python37', 'C:\\Users\\Enzon\\Documents\\Projects\\MEP\\mep-scripts\\venv', 'C:\\Users\\Enzon\\Documents\\Projects\\MEP\\mep-scripts\\venv\\lib\\site-packages', 'C:\\Program Files\\JetBrains\\PyCharm 2020.1.2\\plugins\\python\\helpers\\pycharm_matplotlib_backe#nd', 'C:\\Users\\Enzon\\Documents\\Projects\\MEP\\mep-scripts', 'C:/Users/Enzon/Documents/Projects/MEP/mep-scripts'])

print(sys.path)

#import pip
#help(pip)

#def import_or_install(package):
	#try:
	#	__import__(package)
	#except ImportError:
	#	print('installing '+package)
	#	pip.main(['install', package])

#import_or_install('pip')
#import _imaging
#sys.path.append('C:/Users/Enzon/Documents/Projects/MEP/mep-scripts/venv/Lib/site-packages')

#import_or_install('socket')
#import_or_install('ctypes')
#import_or_install('imageio')
#import_or_install('PIL')
#import_or_install('Pillow')
#import glob
#import imageio

#import PIL.Image
#from PIL import Image
#help(PIL)
#PIL.ImageChops

#sys.path.append('C:\Users\Enzon\Documents\Projects\MEP\mep-scripts\venv\Lib\site-packages')
#import numpy as np
#help(np)



#help(gl)

gl.linewidth(1)
gl.view(64)

# round(0.1)

nRot = 50
# data_path = 'C:/Users/enzon/Documents/Projects/MEP/mep-scripts/Data/Human/Processed'
data_path = os.path.join('/', 'mnt', 'tosh', 'Projects', 'MEP', 'mep-scripts', 'Data', 'Human', 'Processed')
subject_isolation_path_list = glob.glob(os.path.join(data_path, '*', 'Isolations'))
analysis_path = os.path.join('/', 'mnt', 'tosh', 'Projects', 'MEP', 'mep-scripts', 'Data', 'Human', 'Analysis', 'imageSequenceFolders_isolations')
if not os.path.exists(analysis_path):
	os.mkdir(analysis_path)
# analysis_path = 'C:/Users/enzon/Documents/Projects/MEP/mep-scripts/Data/Human/Analysis/imageSequenceFolders'

# subjects = ['patient', 'control1'] #######################################################################

for subject_isolation_path in subject_isolation_path_list:

	subject = subject_isolation_path.split(os.sep)[-2]
	if subject == 'control4':
		## full-ci
		structure_path_list = glob.glob(os.path.join(subject_isolation_path, '*', '*annotation*flirtedRigid*_COMcentered.nii.gz'))
		for structure_path in structure_path_list:
			print(f'structure path = {structure_path}')
			gl.loadimage(structure_path)
			gl.colorname(0, 'FABulous_WT')

			structure_name = structure_path.split(os.sep)[-1].split('_')[-2]
			output_path = os.path.join(analysis_path, subject + '_' + structure_name + '_both')
			if not os.path.exists(output_path):
				os.mkdir(output_path)

			for iRot in range(nRot):
				print(iRot)
				Rot = iRot * (360/nRot)
				print(Rot)
				gl.azimuthelevation(round(Rot), 20)
				filepath = output_path+'/imagenew'+'_'+str(round(Rot)).rjust(3, '0')
				gl.savebmp(filepath)
	elif subject == 'patient':
		## full-ci
		structure_path_list = glob.glob(os.path.join(subject_isolation_path, '*', '*annotation*flirtedRigid*_COMcentered.nii.gz'))
		for structure_path in structure_path_list:
			print(f'structure path = {structure_path}')
			gl.loadimage(structure_path)
			gl.colorname(0, 'FABulous_PAT')

			structure_name = structure_path.split(os.sep)[-1].split('_')[-2]
			output_path = os.path.join(analysis_path, subject + '_' + structure_name + '_both')
			if not os.path.exists(output_path):
				os.mkdir(output_path)

			for iRot in range(nRot):
				print(iRot)
				Rot = iRot * (360/nRot)
				print(Rot)
				gl.azimuthelevation(round(Rot), 20)
				filepath = output_path+'/imagenew'+'_'+str(round(Rot)).rjust(3, '0')
				gl.savebmp(filepath)