#!/bin/bash



# Define function or import packages
normalize () {
  local image=${1}
  local image_normalized=${2}
  local min_max=($(fslstats KO_10/KO_10.nii -R))
  local min=${min_max[0]}
  local max=${min_max[1]}
  local max_adjusted=$(echo "$max-$min" | bc)
  fslmaths $image -sub $min -div $max_adjusted $image_normalized
}



# Define folder and file paths
subject_name="WT_50"
folder_working="/home/enzo/Desktop/Data/Mouse/Processed/${subject_name}" # Alternatively define this by current working directory
folder_atlas_reference="${FSLDIR}/data/standard"
folder_atlas_annotated="${FSLDIR}/data/atlases/AMBMC/"
folder_atlas_xml="${FSLDIR}/data/atlases/"
folder_atlas_LUT="${FSLDIR}/etc/luts/"
folder_atlas_config="${FSLDIR}/etc/flirtsch/"

#~ image="${folder_working}/mouse_mri_data/Pax5_R31Q_mutant/${subject_name}.nii" #~ input?
image="${folder_working}/${subject_name}.nii" #~ input?
image_inmasked="${image%.*}_inmasked.nii"
image_inmasked_flirted="${image_inmasked%.*}_flirted.nii"
image_flirted="${folder_working}/${subject_name}_flirted.nii"
image_inmask="${folder_working}/${subject_name}_mask_t=500_v=380_k=6.mask.nii.gz"
image_inmask_bin="${image_inmask/.*}_bin.mask.nii.gz"
image_inmask_flirted="${folder_working}/${subject_name}_mask_t=500_v=380_k=6_bin_flirted.mask.nii.gz"
image_dir="${image%.*}"
mkdir -p $image_dir
image_name="${image##*/}"
image_name="${image_name%.*}"

image_reference="${folder_atlas_reference}/AMBMC_model.nii.gz" #~ input?
image_reference_name="${image_reference##*/}"
image_reference_name="${image_reference_name%.*}"

image_warpaffine="${image_dir}/${image_name}_to_${image_reference_name}_warpaffine.mat"
image_warpcoef="${image_dir}/${image_name}_to_${image_reference_name}_warpcoef.nii"
image_warped="${image_dir}/${image_name}_to_${image_reference_name}_warped.nii"
image_warpcoef_inverted="${image_dir}/${image_reference_name}_to_${image_name}_warpcoef_inverted.nii"

image_basalganglia="${folder_atlas_annotated}/AMBMC-c57bl6-basalganglia-labels-15um.nii.gz"
image_cerebellum="${folder_atlas_annotated}/AMBMC-c57bl6-cerebellum-labels-15um.nii.gz"
image_cortex="${folder_atlas_annotated}/AMBMC-c57bl6-cortex-labels-15um.nii.gz"
image_hippocampus="${folder_atlas_annotated}/AMBMC-c57bl6-hippocampus-labels-15um.nii.gz"
image_structures=($image_basalganglia $image_cerebellum $image_cortex $image_hippocampus)
#~ image_structures=($image_cerebellum)

image_basalganglia_invwarped="${image_dir}/basalganglia_to_${image_name}_invwarped.nii"
image_cerebellum_invwarped="${image_dir}/cerebellum_to_${image_name}_invwarped.nii"
image_cortex_invwarped="${image_dir}/cortex_to_${image_name}_invwarped.nii"
image_hippocampus_invwarped="${image_dir}/hippocampus_to_${image_name}_invwarped.nii"
image_structures_invwarped=($image_basalganglia_invwarped $image_cerebellum_invwarped $image_cortex_invwarped $image_hippocampus_invwarped)
#~ image_structures_invwarped=($image_cerebellum_invwarped)



echo "Working on ${subject_name}"

echo .
echo "normalize ${image_inmasl}"
normalize $image_inmask $image_inmask_bin

echo mul
fslmaths $image -mul $image_inmask_bin $image_inmasked

echo .
echo flirt
flirt 	-in $image_inmasked \
		-ref $image_reference \
		-omat $image_warpaffine \
		-out $image_inmasked_flirted \
		-verbose 1
		
echo "flirt actual image"
flirt 	-in $image \
		-ref $image_reference \
		-out $image_flirted \
		-init $image_warpaffine \
		-applyxfm \
		-verbose 1

#~ echo .
#~ echo fnirt
#~ fnirt 	--config=AMBMC_config.cnf \
		#~ --in=$image \
		#~ --ref=$image_reference \
		#~ --cout=$image_warpcoef \
		#~ --aff=$image_warpaffine \
		#~ --inmask=$image_inmask \
		#~ --applyinmask=1 \
		#~ --iout=$image_warped \
		#~ --verbose
		
		
		
#~ --miter=1,1,1,1,1,1 \

#~ applywarp 	--in=$image \
			#~ --ref=$image_reference \
			#~ --warp=$image_warpcoef \
			#~ --out=$image_warped \
			#~ --verbose
			
			

#~ echo .
#~ echo invwarp
#~ invwarp --ref=$image \
		#~ --warp=$image_warpcoef \
		#~ --out=$image_warpcoef_inverted \
		#~ --verbose
		
#~ echo .
#~ echo "inverted applywarps"
#~ for i in ${!image_structures[@]}; do
	#~ applywarp 	--ref=$image \
				#~ --in=${image_structures[i]} \
				#~ --warp=$image_warpcoef_inverted \
				#~ --out=${image_structures_invwarped[i]} \
				#~ --interp=nn \
				#~ --verbose
#~ done
