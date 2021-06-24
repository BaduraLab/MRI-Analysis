#!/bin/bash

for fp_in in $(find /home/enzo/Desktop/Data/Mouse/Raw/*nii); do
	fp_out="/home/enzo/Desktop/Data/Mouse/Processed/${fp_in##*/}"
	fp_out="${fp_out%.*}.nrrd"
	fp_out_nii="${fp_out%.*}.nii"
	python /home/enzo/Desktop/mri_convert_itk.py $fp_in $fp_out
	python /home/enzo/Desktop/mri_convert_itk.py $fp_out $fp_out_nii
done
	
