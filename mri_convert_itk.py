#!/usr/bin/env python

import itk, sys, os
    
if __name__ == '__main__':
	# Parse arguments
	image_input  = sys.argv[1]
	image_output = sys.argv[2]
	
	image = itk.imread(image_input)
	
	itk.imwrite(image, image_output)
