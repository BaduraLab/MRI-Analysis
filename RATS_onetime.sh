#!/bin/bash

# Define parameters
folder_working="/home/enzo/Desktop/WT_50_RATS_MM_test"
folder_RATS="/home/enzo/Desktop/rats/distribution"
folder_slicer="/usr/local/slicer"
image="$folder_working/WT_50.nii" # (rat_t1.nrrd) Input image
image_name="${image%.*}"
image_nrrd="${image_name}.nrrd"
image_mask="/home/enzo/Desktop/C2-Composite_bin.nii"
image_vtp="/home/enzo/Desktop/C2-Composite_bin_vtp.vtp"
image_vtp_mask="/home/enzo/Desktop/C2-Composite_bin_vtp_mask.nrrd"
a=0
n=1000
image_vtp_mask_nii="/home/enzo/Desktop/C2-Composite_bin_vtp_mask_a=${a}_n=${n}.nii"

echo RATS_LOGISMOS
${folder_RATS}/RATS_LOGISMOS 	$image_nrrd \
								$image_mask \
								$image_vtp \
								-a ${a} -n ${n}
								
echo MeshToLabelMap
${folder_slicer}/Slicer --launch MeshToLabelMap \
	--reference_volume $image_mask \
	--input_mesh $image_vtp \
	--output_labelmap $image_vtp_mask
	
# Convert .nrrd output to .nii with custom slicer python script
echo "Convert nrrd mask to nii mask"
${folder_slicer}/Slicer --no-main-window --no-splash --python-script ${folder_working}/mri_convert_slicer.py $image_vtp_mask $image_vtp_mask_nii
