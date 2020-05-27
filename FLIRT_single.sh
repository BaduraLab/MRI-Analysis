#!/bin/bash

#~ Example:
#~ script -c "bash -x Desktop/mep-scripts/FLIRT_single.sh WT_50"

start=`date +%h`



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
both(){ ( echo "${@:2}" && "${@:2}" ) | tee "$1" ;}



## Define folder and file paths
#~ subject_name="WT_50"
subject_name=${1}

folder_working="/home/enzo/Desktop/Data/Mouse/Processed_New/${subject_name}" # Alternatively define this by current working directory
folder_FLIRT="${folder_working}/FLIRT"
mkdir -p $folder_FLIRT
cp /home/enzo/Desktop/mep-scripts/FLIRT_single.sh $folder_FLIRT/

# Log files
FLIRT_log="${folder_FLIRT}/FLIRT_log"

# Folders of reference
folder_atlas_reference="${FSLDIR}/data/standard/allen_new"
folder_atlas_annotated="${FSLDIR}/data/atlases/AMBMC"
folder_atlas_xml="${FSLDIR}/data/atlases"
folder_atlas_LUT="${FSLDIR}/etc/luts"
folder_atlas_config="${FSLDIR}/etc/flirtsch"

# Define image path and reorient, use reoriented image as image
image="${folder_working}/${subject_name}.nii" # input?
image_reoriented="${folder_FLIRT}/${subject_name}.nii"
fslreorient2std $image $image_reoriented
image=$image_reoriented

# Define image mask path and reorient, use reoriented mask as mask
image_inmask="${folder_working}/${subject_name}_mask_t=500_v=380_k=6.mask.nii.gz"
image_inmask_reoriented="${folder_FLIRT}/${subject_name}_mask_t=500_v=380_k=6.mask.nii.gz"
fslreorient2std $image_inmask $image_inmask_reoriented
image_inmask=$image_inmask_reoriented

# Define image and mask FLIRT output paths
image_inmasked="${image%.*}_inmasked.nii"
image_inmasked_flirted="${image_inmasked%.*}_flirted.nii"
image_flirted="${folder_FLIRT}/${subject_name}_flirted.nii"
image_inmask_bin="${image_inmask/.*}_bin.mask.nii.gz"
image_inmask_flirted="${folder_FLIRT}/${subject_name}_mask_t=500_v=380_k=6_bin_flirted.mask.nii.gz"
image_name="${image##*/}"
image_name="${image_name%.*}"

# Define image reference path and check whether it is reoriented to reference
image_reference="${folder_atlas_reference}/average_template_25_to_AMBMC_flirted.nii.gz" # input?
image_reference_name="allen_model" # input?
echo "Print reorientation-of-reference-to-standard-reference matrix"
fslreorient2std $image_reference # Prints matrix, if identity matrix then ok

# Rescale image to have same value range as reference, assuming minimum of reference is 0
image_reference_min_max=($(fslstats $image_reference -R))
image_reference_max=${image_reference_min_max[1]}
normalize $image $image
fslmaths $image -mul $image_reference_max $image

# Define FLIRT paths
image_warpaffine="${folder_FLIRT}/${image_name}_to_${image_reference_name}_warpaffine.mat"
image_warpaffine_inverted="${folder_FLIRT}/${image_reference_name}_to_${image_name}_warpaffine_inverted.mat"

# Define annotation volume paths
image_allen_annotation="${folder_atlas_reference}/annotation_25_to_AMBMC_flirted.nii.gz" # input?
image_allen_annotation_invflirted="${folder_FLIRT}/allen_annotation_to_${image_name}_invwarped.nii"



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
both $FLIRT_log \
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
		


# Invert FLIRT (reference2native)
echo
echo "Invert FLIRT"
convert_xfm -omat $image_warpaffine_inverted -inverse $image_warpaffine

# FLIRT applywarp (reference2native)
echo
echo "Invert FLIRT allen to image"
flirt 	-in $image_allen_annotation \
		-ref $image \
		-out $image_allen_annotation_invflirted \
		-init $image_warpaffine_inverted \
		-applyxfm \
		-interp nearestneighbour \
		-verbose 1
		
		
		
end=`date +%h`
runtime=$((end-start))
echo $runtime
