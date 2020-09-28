import gl
import sys
print(sys.version)
print(gl.version())
import os

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

round(0.1)
nRot = 50
data_path = 'C:/Users/enzon/Documents/Projects/MEP/mep-scripts/Data/Human/Processed'
analysis_path = 'C:/Users/enzon/Documents/Projects/MEP/mep-scripts/Data/Human/Analysis/imageSequenceFolders'

subjects = ['patient', 'control1'] #######################################################################

for subject in subjects:

	for i in range(4):

		if i==0:
			## full-ci
			input_path = data_path + '/' + subject+'/'+subject+'_reoriented.nii.gz'
			overlay_path = data_path + '/' + subject+'/'+subject+'_annotation_orsuit_thrarg_adjusted_lobular_mc_ci.nii.gz'
			output_path = analysis_path + '/' + subject+'_50rot_overlay_full_ci'
		elif i==1:
			## ci-ci
			input_path = data_path + '/' + subject+'/'+subject+'_reoriented_ci.nii.gz'
			overlay_path = data_path + '/' + subject+'/'+subject+'_annotation_orsuit_thrarg_adjusted_lobular_mc_ci.nii.gz'
			output_path = analysis_path + '/' + subject+'_50rot_overlay_ci_ci'
		elif i==2:
			## full-si, only cutout interesting really
			input_path = data_path + '/' + subject+'/'+subject+'_reoriented.nii.gz'
			overlay_path = data_path + '/' + subject+'/'+subject+'_annotation_orsuit_thrarg_adjusted_lobular_mc_ci.nii.gz'
			overlay2_path = data_path + '/' + subject+'/'+subject+'_annotation_subcortical_thrarg_si.nii.gz'
			output_path = analysis_path + '/' + subject+'_50rot_overlay_full_si'
		else:
			## si-si
			input_path = data_path + '/' + subject+'/'+subject+'_reoriented_si.nii.gz'
			overlay_path = data_path + '/' + subject+'/'+subject+'_annotation_subcortical_thrarg_si.nii.gz'
			output_path = analysis_path + '/' + subject+'_50rot_overlay_si_si'

		print(input_path)
		print(overlay_path)

		for j in range(2):

			if j==1:
				output_path = output_path+'_cutout'
			output_gif_path = output_path+'.gif'
			print(output_path)

			if not os.path.exists(output_path):
				os.mkdir(output_path)

			if i==2:
				gl.loadimage(input_path)
				gl.overlayload(overlay_path)
				gl.overlayload(overlay2_path)
				gl.opacity(1, 40)
				gl.opacity(2, 70)
				gl.minmax(1, 0, 34)
				gl.minmax(2, 6, 10)
				gl.colorname(1, 'x_rain')
				gl.colorname(2, 'x_rain')
				#gl.fullscreen(0)
			else:
				gl.loadimage(input_path)
				gl.overlayload(overlay_path)
				gl.opacity(1, 40)
				gl.minmax(1, 0, 34)
				gl.colorname(1, 'x_rain')
				#gl.fullscreen(0)

			if j==1:
				if i<2:
					gl.cutout(0.5, 0.5, 0.5, 0, 0, 1)
				else:
					gl.cutout(0.53, 0.42, 0.39, 0, 0, 1)

			# help(gl.colorname)

			gl.shadername('overlaysurface')

			#im_list = []
			for iRot in range(nRot):
				print(iRot)
				Rot = iRot * (360/nRot)
				print(Rot)
				gl.azimuthelevation(round(Rot), 20)
				filepath = output_path+'/imagenew'+'_'+str(round(Rot)).rjust(3, '0')
				gl.savebmp(filepath)
				#im = PIL.Image.open(filepath)
				#im_list.append(im)
				
			#print('saving '+output_gif_path)
			#im_list[0].save(output_gif_path, save_all=True, append_images=images[1:], optimize=False, duration=40, loop=0)

			#files = glob.glob(os.path.join(output_path, '*'))
			#os.listdir(output_path)
			#files = [f"{folder}/{file}" for file in os.listdir(output_path)]
			##images = [imageio.imread(file) for file in files]
			
			##imageio.mimwrite(os.path.join(analysis_path, 'movie.gif'), images, fps=20)






	# gl.resetdefaults()
	# gl.viewcoronal(True)
	#gl.view(2)
	#gl.orthoviewmm(0.5,0.5,0.5)