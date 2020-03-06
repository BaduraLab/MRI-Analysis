#!/bin/bash

find . -type f -name "*.nii*"

for imagename in $(find . -type f -name "*.nii*"); do

	echo "Start image reorientation for ${imagename}"

	fslorient -setqformcode 1 $imagename
	fslorient -setqform 0.05 0 0 0 0 -0.05 0 0 0 0 0.05 0 0 0 0 1 $imagename
	fslorient -setsformcode 0 $imagename
	fslreorient2std $imagename

	fslhd $imagename
	#~ fsleyes $imagename

done
