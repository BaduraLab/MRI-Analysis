#!/bin/bash



# Define function or import packages
normalize () {
  local image=${1}
  local image_normalized=${2}
  local min_max=($(fslstats $image -R))
  local min=${min_max[0]}
  local max=${min_max[1]}
  local max_adjusted=$(echo "$max-$min" | bc)
  fslmaths $image -sub $min -div $max_adjusted $image_normalized
}



## Define folder and file paths
#~ subject_name="WT_50"
subject_name=${1}
folder_working="/home/enzo/Desktop/Data/Mouse/Processed/${subject_name}" # Alternatively define this by current working directory

# Folders of reference
folder_atlas_reference="${FSLDIR}/data/standard"
folder_atlas_annotated="${FSLDIR}/data/atlases/AMBMC/"
folder_atlas_xml="${FSLDIR}/data/atlases/"
folder_atlas_LUT="${FSLDIR}/etc/luts/"
folder_atlas_config="${FSLDIR}/etc/flirtsch/"

# Define image path and reorient, use reoriented image as image
image="${folder_working}/${subject_name}.nii" #~ input?
image_reoriented="${image/.*}_reoriented"
fslreorient2std $image $image_reoriented
image=$image_reoriented

# Define image mask path and reorient, use reoriented mask as mask
image_inmask="${folder_working}/${subject_name}_mask_t=500_v=380_k=6.mask.nii.gz"
image_inmask_reoriented="${image_inmask/.*}_reoriented.nii"
fslreorient2std $image_inmask $image_inmask_reoriented
image_inmask=$image_inmask_reoriented

# Define image and mask FLIRT output paths
image_inmasked="${image%.*}_inmasked.nii"
image_inmasked_flirted="${image_inmasked%.*}_flirted.nii"
image_flirted="${folder_working}/${subject_name}_flirted.nii"
image_inmask_bin="${image_inmask/.*}_bin.mask.nii.gz"
image_inmask_flirted="${folder_working}/${subject_name}_mask_t=500_v=380_k=6_bin_flirted.mask.nii.gz"
image_name="${image##*/}"
image_name="${image_name%.*}"

# Define image reference path and check whether it is reoriented to reference
image_reference="${folder_atlas_reference}/AMBMC_model_reoriented.nii.gz" #~ input?
image_reference_name="${image_reference##*/}"
image_reference_name="${image_reference_name%%.*}"
echo "Print reorientation-of-reference-to-standard-reference matrix"
fslreorient2std $image_reference # Prints matrix, if identity then ok

# Define FLIRT and FNIRT paths
image_warpaffine="${folder_working}/${image_name}_to_${image_reference_name}_warpaffine.mat"
image_warpcoef="${folder_working}/${image_name}_to_${image_reference_name}_warpcoef.nii"
image_warped="${folder_working}/${image_name}_to_${image_reference_name}_warped.nii"
image_warpaffine_inverted="${folder_working}/${image_reference_name}_to_${image_name}_warpaffine_inverted.mat"
image_warpcoef_inverted="${folder_working}/${image_reference_name}_to_${image_name}_warpcoef_inverted.nii"

# Define labelled reference input paths
image_basalganglia="${folder_atlas_annotated}/AMBMC-c57bl6-basalganglia-labels-15um_reoriented.nii.gz"
image_cerebellum="${folder_atlas_annotated}/AMBMC-c57bl6-cerebellum-labels-15um_reoriented.nii.gz"
image_cortex="${folder_atlas_annotated}/AMBMC-c57bl6-cortex-labels-15um_reoriented.nii.gz"
image_hippocampus="${folder_atlas_annotated}/AMBMC-c57bl6-hippocampus-labels-15um_reoriented.nii.gz"
image_structures=($image_reference $image_basalganglia $image_cerebellum $image_cortex $image_hippocampus)

# Define labelled reference outputs paths which are warped to native space
image_reference_invwarped="${folder_working}/model_to_${image_name}_invwarped.nii.gz"
image_basalganglia_invwarped="${folder_working}/basalganglia_to_${image_name}_invwarped.nii"
image_cerebellum_invwarped="${folder_working}/cerebellum_to_${image_name}_invwarped.nii"
image_cortex_invwarped="${folder_working}/cortex_to_${image_name}_invwarped.nii"
image_hippocampus_invwarped="${folder_working}/hippocampus_to_${image_name}_invwarped.nii"
image_structures_invwarped=($image_reference_invwarped $image_basalganglia_invwarped $image_cerebellum_invwarped $image_cortex_invwarped $image_hippocampus_invwarped)



# FLIRT
echo "Working on ${subject_name}"

echo
echo "normalize ${image_inmask}"
normalize $image_inmask $image_inmask_bin

echo
echo mul
fslmaths $image -mul $image_inmask_bin $image_inmasked



# FLIRT (native2reference)
echo
echo flirt
flirt 	-in $image_inmasked \
		-ref $image_reference \
		-omat $image_warpaffine \
		-out $image_inmasked_flirted \
		-verbose 1
		
# FLIRT applywarp to image without mask (native2reference)
echo
echo "flirt actual image"
flirt 	-in $image \
		-ref $image_reference \
		-out $image_flirted \
		-init $image_warpaffine \
		-applyxfm \
		-verbose 1



# FNIRT (native2reference)
echo
echo fnirt
fnirt 	--config=AMBMC_config.cnf \
		--in=$image \
		--ref=$image_reference \
		--cout=$image_warpcoef \
		--aff=$image_warpaffine \
		--inmask=$image_inmask_bin \
		--applyinmask=1 \
		--iout=$image_warped \
		--verbose
	


# Invert FNIRT (reference2native) $image_flirted could be the one
echo
echo "Invert FNIRT"
invwarp --ref=$image_flirted \
		--warp=$image_warpcoef \
		--out=$image_warpcoef_inverted \
		--verbose
	
		

# Warp labelled references to native space
echo
echo "inverted applywarps"
for i in ${!image_structures[@]}; do
	applywarp 	--ref=$image \
				--in=${image_structures[i]} \
				--warp=$image_warpcoef_inverted \
				--out=${image_structures_invwarped[i]} \
				--interp=nn \
				--verbose
done



# FNIRT reference directly to native space



#~ --postmat=$image_warpaffine_inverted \

#~ --miter=1,1,1,1,1,1 \

#~ applywarp 	--in=$image \
			#~ --ref=$image_reference \
			#~ --warp=$image_warpcoef \
			#~ --out=$image_warped \
			#~ --verbose
			
#~ # Invert FLIRT (reference2native)
#~ echo
#~ echo "Invert FLIRT"
#~ convert_xfm -omat $image_warpaffine_inverted -inverse $image_warpaffine

#~ image_structures=($image_cerebellum)
#~ image_structures_invwarped=($image_cerebellum_invwarped)

#~ image_dir="${image%.*}"
#~ mkdir -p $image_dir
