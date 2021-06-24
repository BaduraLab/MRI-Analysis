import vtk, sys, os

def readnifti(filename):
    """Read image in nrrd format."""
    reader = vtk.vtkNIFTIImageReader()
    reader.SetFileName(filename)
    reader.Update()
    info = reader.GetInformation()
    return reader.GetOutput(), info
    
def readnrrd(filename):
    """Read image in nrrd format."""
    reader = vtk.vtkNrrdReader()
    reader.SetFileName(filename)
    reader.Update()
    info = reader.GetInformation()
    return reader.GetOutput(), info

def writenifti(image,filename, info):
    """Write nifti file."""
    writer = vtk.vtkNIFTIImageWriter()
    writer.SetInputData(image)
    writer.SetFileName(filename)
    writer.SetInformation(info)
    writer.Write()
    
def writenrrd(image,filename, info):
    """Write nifti file."""
    writer = vtk.vtkNrrdWriter()
    writer.SetInputData(image)
    writer.SetFileName(filename)
    writer.SetInformation(info)
    writer.Write()
    
if __name__ == '__main__':
	# Parse arguments
	image_input  = sys.argv[1]
	image_output = sys.argv[2]
	
	# Get input image file extension
	_, input_extension = os.path.splitext(image_input)

	print(input_extension)

	# Read and write data
	if input_extension == '.nii':
		m, info = readnifti(image_input)
		writenrrd(m, image_output, info)
		print(image_output)
	elif input_extension == '.nrrd':
		m, info = readnrrd(image_input)
		writenifti(m, image_output, info)
		print(image_output)
