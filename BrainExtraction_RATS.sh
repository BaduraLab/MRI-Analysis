#!/bin/bash
# Use RATS to generate brain mask images for all files within input folder
#
# Usage: BrainExtraction_RATS.sh "data folder"
#
# Use Rapid Automatic Tissue Segmentation to generate a VTK surface model (.vtp file)
# Convert that surface model into a brain mask using 3D slicer.
#
# Requires the installation of RATS and 3D slicer


start=`date +%s`



# Define parameters
folder_data=${1}
# folder_RATS=${2}
# folder_RATS="/home/enzo/Desktop/BrainExtraction/rats/distribution" # Could be input (RATS installation folder)
# folder_slicer="/usr/local/slicer"

for image_dir in $(find $folder_data/* -prune -type d); do
	echo "Working on ${image_dir}"
	
	image="${image_dir}/${image_dir##*/}.nii.gz"
	image_name="${image%%.*}"
	#~ image_nrrd="${image_name}.nrrd"
	t=(500) # Intensity threshold [500]
	v=(380) # Volume threshold [380]
	k=(6) # Diameter of structuring element size K1 [6]
	a=(0) # Weight factor between the two cost terms [0]
	n=(2000) # Number of vertices which describe the surface [2000]

	# Multiply image by 208500 which is approximately the multiplication factor between Bru2Nii output and previous imagej plugin output
	fslmaths $image -mul 208500 $image

	echo "enter parameter for loop"
	for i in ${!t[@]}; do
		echo "t=${t[i]}"
		echo "v=${v[i]}"
		echo "k=${k[i]}"
		echo "a=${a[i]}"
		echo "n=${n[i]}"

		image_mask="${image_name}_mask_t=${t[i]}_v=${v[i]}_k=${k[i]}.nii.gz" # (rat_t1_mm.nii.gz) Pre-segmentation image
		image_vtp="${image_name}_vtp_t=${t[i]}_v=${v[i]}_k=${k[i]}_a=${a[i]}_n=${n[i]}.vtp" # (rat_t1_logismos.vtp) Output mesh (.vtp)
		image_vtp_mask="${image_name}_vtp_mask_t=${t[i]}_v=${v[i]}_k=${k[i]}_a=${a[i]}_n=${n[i]}.nii.gz" # Final mask
		
		# RATS_MM
		if test -f "$image_mask"; then
			echo "RATS_MM unneccessary as output file already exists"
		else			
			echo RATS_MM
			#${folder_RATS}/RATS_MM 	$image \
			RATS_MM 	$image \
								$image_mask \
								-t ${t[i]} -v ${v[i]} -k ${k[i]}
		fi
		
		# RATS_LOGISMOS
		if test -f "$image_vtp"; then
			echo "RATS_LOGISMOS unneccessary as output file already exists"
		else	
			echo RATS_LOGISMOS
			#${folder_RATS}/RATS_LOGISMOS 	$image \
			RATS_LOGISMOS 	$image \
											$image_mask \
											$image_vtp \
											-a ${a[i]} -n ${n[i]}
		fi	
							
		# MeshToLabelMap
		if test -f "$image_vtp_mask"; then
			echo "MeshToLabelMap unneccessary as output file already exists"
		else		
			echo MeshToLabelMap
			#${folder_slicer}/Slicer --launch MeshToLabelMap \
			Slicer --launch MeshToLabelMap \
						 --reference_volume $image_mask \
						 --input_mesh $image_vtp \
						 --output_labelmap $image_vtp_mask
		fi
		
	done
	
done



end=`date +%s`
runtime=$((end-start))
echo $runtime
