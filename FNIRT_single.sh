#!/bin/bash

#~ Example:
#~ script -c "bash -x Desktop/mep-scripts/FNIRT_single.sh WT_50"

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
datestr="$(date +%F)-$(date +%H)-$(date +%M)-$(date +%S)"
folder_FNIRT="${folder_working}/FNIRT_${datestr}"
mkdir -p $folder_FNIRT
cp /home/enzo/Desktop/mep-scripts/FNIRT_single.sh $folder_FNIRT/

# Log files
FNIRT_log="${folder_FNIRT}/FNIRT_log"

# Folders of reference
folder_atlas_reference="${FSLDIR}/data/standard/allen/"
folder_atlas_annotated="${FSLDIR}/data/atlases/AMBMC/"
folder_atlas_xml="${FSLDIR}/data/atlases/"
folder_atlas_LUT="${FSLDIR}/etc/luts/"
folder_atlas_config="${FSLDIR}/etc/flirtsch/"

# Define image path and reorient, use reoriented image as image
image="${folder_FLIRT}/${subject_name}.nii" # input?

# Define image mask path and reorient, use reoriented mask as mask
image_inmask="${folder_FLIRT}/${subject_name}_mask_t=500_v=380_k=6.mask.nii.gz"

# Define image and mask FLIRT output paths
image_inmask_bin="${image_inmask/.*}_bin.mask.nii.gz"
image_inmask_bin_warped="${image_inmask/.*}_bin.mask_warped.nii.gz"
image_inmask_bin_warpcoef="${image_inmask/.*}_bin.mask_warpcoef.nii.gz"
image_name="${image##*/}"
image_name="${image_name%.*}"

# Define image reference path and check whether it is reoriented to reference
image_reference="${folder_atlas_reference}/average_template_25_to_AMBMC_flirted.nii.gz" # input?
image_reference_name="allen_model" # input?
echo "Print reorientation-of-reference-to-standard-reference matrix"
fslreorient2std $image_reference # Prints matrix, if identity then ok

# Define FLIRT and FNIRT paths
image_warpaffine="${folder_FLIRT}/${image_name}_to_${image_reference_name}_warpaffine.mat"
image_warpcoef="${folder_FNIRT}/${image_name}_to_${image_reference_name}_warpcoef.nii"
image_warped="${folder_FNIRT}/${image_name}_to_${image_reference_name}_warped.nii"
image_warpaffine_inverted="${folder_FLIRT}/${image_reference_name}_to_${image_name}_warpaffine_inverted.mat"
image_warpcoef_inverted="${folder_FNIRT}/${image_reference_name}_to_${image_name}_warpcoef_inverted.nii"

# Define annotation volume paths
image_allen_annotation="${folder_atlas_reference}/annotation_25_to_AMBMC_flirted.nii.gz" # input?
image_allen_annotation_bin="${folder_atlas_reference}/annotation_25_reoriented_bin_to_AMBMC_flirted.nii.gz"
image_allen_annotation_invwarped="${folder_FNIRT}/allen_annotation_to_${image_name}_invwarped.nii" # input?
image_allen_annotation_invflirted="${folder_FLIRT}/allen_annotation_to_${image_name}_invflirted.nii"



# FNIRT initial guess by warping mask to reference(native2reference)echo
echo fnirt
both $FNIRT_log \
fnirt 	--config=AMBMC_config.cnf \
		--in=$image_inmask_bin \
		--ref=$image_reference \
		--cout=$image_inmask_bin_warpcoef \
		--aff=$image_warpaffine \
		--iout=$image_inmask_bin_warped \
		--verbose \
		--subsamp=4,4 \
		--miter=5,5 \
		--infwhm=0.18,0.15 \
		--reffwhm=0.06,0.06 \
		--lambda=150,50 \
		--estint=0,0 \
		--applyrefmask=0,0 \
		--ssqlambda=0

# FNIRT (native2reference)
echo
echo fnirt
both $FNIRT_log \
fnirt 	--config=AMBMC_config.cnf \
		--in=$image \
		--ref=$image_reference \
		--cout=$image_warpcoef \
		--inwarp=$image_inmask_bin_warpcoef \
		--inmask=$image_inmask_bin \
		--applyinmask=1,1,1,1,1,1 \
		--iout=$image_warped \
		--verbose \
		--subsamp=4,4,2,2,1,1 \
		--miter=5,5,5,5,5,10 \
		--infwhm=0.18,0.15,0.12,0.09,0.06,0 \
		--reffwhm=0.06,0.06,0.03,0.03,0,0 \
		--lambda=150,50,30,15,5,1 \
		--estint=1,1,1,1,1,0 \
		--applyrefmask=0,0,0,0,0,0 \
		--ssqlambda=0
	
		#~ --aff=$image_warpaffine \	
#~ --infwhm=1.44,0.72\
#~ --reffwhm=0.96,0.48\
#~ --lambda=75,5 \
#~ --estint=1 \
#~ --applyrefmask=0 \
#~ --ssqlambda=0

#~ --subsamp=4,4,2,2,1,1
#~ --miter=5,5,5,5,5,10
#~ --infwhm=0.18,0.15,0.12,0.09,0.06,0
#~ --reffwhm=0.06,0.06,0.03,0.03,0,0
#~ --lambda=150,75,50,30,20,15
#~ --estint=1,1,1,1,1,0
#~ --applyrefmask=0,0,0,0,0,0

#~ --warpres=0.66,1,0.46
#~ --biasres=3,5,2



# Invert FNIRT (reference2native)
echo
echo "Invert FNIRT"
invwarp --ref=$image_flirted \
		--warp=$image_warpcoef \
		--out=$image_warpcoef_inverted \
		--verbose

# FNIRT applywarp (reference2native)
# Warp labelled references to native space
echo
echo "inverted applywarp"
applywarp 	--ref=$image \
			--in=$image_allen_annotation \
			--warp=$image_warpcoef_inverted \
			--out=$image_allen_annotation_invwarped \
			--interp=nn \
			--verbose



end=`date +%h`
runtime=$((end-start))
echo $runtime
