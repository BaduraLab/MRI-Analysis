#!/bin/bash

find . -type f -name "*.nii*"

for imagename in $(find . -type f -name "*.nii*"); do

	echo "Start image reorientation for ${imagename}"
	
	imagename=${imagename##*/}
	imagename_reoriented="${imagename%%.*}_reoriented.nii.gz"

	sudo cp $imagename $imagename_reoriented
	sudo chmod 777 $imagename_reoriented

	
	
	echo
	echo "header before"
	fslhd $imagename
	
	fslorient -setqformcode 1 $imagename_reoriented
	fslorient -setqform ${1} 0 0 0 0 -${1} 0 0 0 0 ${1} 0 0 0 0 1 $imagename_reoriented
	fslorient -setsformcode 0 $imagename_reoriented
	
	echo
	echo "header after"
	fslhd $imagename_reoriented
	
	
	
	echo
	echo "reorientation matrix to MNI before"
	fslreorient2std $imagename_reoriented
	
	fslreorient2std $imagename_reoriented $imagename_reoriented
	
	echo
	echo "reorientation matrix to MNI after"
	fslreorient2std $imagename_reoriented
	
	

	fsleyes $imagename_reoriented

done
