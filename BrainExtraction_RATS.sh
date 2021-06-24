#!/bin/bash



start=`date +%s`



# Define parameters
folder_working="/home/enzo/Desktop/BrainExtraction"
#~ folder_data="/home/enzo/Desktop/Data/Mouse/Processed" # Could be input
folder_data=${1} # Could be input
folder_RATS="${folder_working}/rats/distribution"
folder_slicer="/usr/local/slicer"

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

	#~ # Convert .nii.gz input to .nrrd with custom slicer python script if nrrd file does not already exists
	#~ if test -f "$image_nrrd"; then
		#~ echo "nrrd already exists, no need to convert nii.gz"
	#~ else
		#~ echo "convert nii.gz to nrrd"
		#~ python Desktop/mep-scripts/mri_convert_itk.py $image $image_nrrd
	#~ fi

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
			${folder_RATS}/RATS_MM 	$image \
									$image_mask \
									-t ${t[i]} -v ${v[i]} -k ${k[i]}
		fi
		
		# RATS_LOGISMOS
		if test -f "$image_vtp"; then
			echo "RATS_LOGISMOS unneccessary as output file already exists"
		else	
			echo RATS_LOGISMOS
			${folder_RATS}/RATS_LOGISMOS 	$image \
											$image_mask \
											$image_vtp \
											-a ${a[i]} -n ${n[i]}
		fi	
							
		# MeshToLabelMap
		if test -f "$image_vtp_mask"; then
			echo "MeshToLabelMap unneccessary as output file already exists"
		else		
			echo MeshToLabelMap
			${folder_slicer}/Slicer --launch MeshToLabelMap \
									--reference_volume $image_mask \
									--input_mesh $image_vtp \
									--output_labelmap $image_vtp_mask
		fi
		
	done
	
done



end=`date +%s`
runtime=$((end-start))
echo $runtime
