#!/bin/bash

#~ Example:
#~ script -c "bash -x Desktop/mep-scripts/FLIRT_and_FNIRT_single.sh WT_50"

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
datestr="$(date +%F)-$(date +%H)-$(date +%M)-$(date +%S)"
folder_FNIRT="${folder_working}/FNIRT_${datestr}"
mkdir -p $folder_FNIRT
cp /home/enzo/Desktop/mep-scripts/FNIRT_single.sh $folder_FNIRT/

# Log files
FLIRT_log="${folder_FNIRT}/FLIRT_log"
FNIRT_log="${folder_FNIRT}/FNIRT_log"

#~ exec 1 2>&1 | tee ${LOG_FILE}
#~ exec > "${folder_FNIRT}/log.txt" 2>&1
#~ set -x

# Folders of reference
folder_atlas_reference="${FSLDIR}/data/standard/allen/"
folder_atlas_annotated="${FSLDIR}/data/atlases/AMBMC/"
folder_atlas_xml="${FSLDIR}/data/atlases/"
folder_atlas_LUT="${FSLDIR}/etc/luts/"
folder_atlas_config="${FSLDIR}/etc/flirtsch/"

# Define image path and reorient, use reoriented image as image
image="${folder_working}/${subject_name}.nii" # input?
image_reoriented="${folder_FNIRT}/${subject_name}.nii"
fslreorient2std $image $image_reoriented
image=$image_reoriented

# Define image mask path and reorient, use reoriented mask as mask
image_inmask="${folder_working}/${subject_name}_mask_t=500_v=380_k=6.mask.nii.gz"
image_inmask_reoriented="${folder_FNIRT}/${subject_name}_mask_t=500_v=380_k=6.mask.nii.gz"
fslreorient2std $image_inmask $image_inmask_reoriented
image_inmask=$image_inmask_reoriented

# Define image and mask FLIRT output paths
image_inmasked="${image%.*}_inmasked.nii"
image_inmasked_flirted="${image_inmasked%.*}_flirted.nii"
image_flirted="${folder_FNIRT}/${subject_name}_flirted.nii"
image_inmask_bin="${image_inmask/.*}_bin.mask.nii.gz"
image_inmask_bin_warped="${image_inmask/.*}_bin.mask_warped.nii.gz"
image_inmask_bin_warpcoef="${image_inmask/.*}_bin.mask_warpcoef.nii.gz"
image_inmask_flirted="${folder_FNIRT}/${subject_name}_mask_t=500_v=380_k=6_bin_flirted.mask.nii.gz"
image_name="${image##*/}"
image_name="${image_name%.*}"

# Define image reference path and check whether it is reoriented to reference
image_reference="${folder_atlas_reference}/average_template_25_to_AMBMC_flirted.nii.gz" # input?
image_reference_name="allen_model" # input?
echo "Print reorientation-of-reference-to-standard-reference matrix"
fslreorient2std $image_reference # Prints matrix, if identity then ok

# Rescale image to have same value range as reference, assuming minimum of reference is 0
image_reference_min_max=($(fslstats $image_reference -R))
image_reference_max=${image_reference_min_max[1]}
normalize $image $image
fslmaths $image -mul $image_reference_max $image

# Define FLIRT and FNIRT paths
image_warpaffine="${folder_FNIRT}/${image_name}_to_${image_reference_name}_warpaffine.mat"
image_warpcoef="${folder_FNIRT}/${image_name}_to_${image_reference_name}_warpcoef.nii"
image_warped="${folder_FNIRT}/${image_name}_to_${image_reference_name}_warped.nii"
image_warpaffine_inverted="${folder_FNIRT}/${image_reference_name}_to_${image_name}_warpaffine_inverted.mat"
image_warpcoef_inverted="${folder_FNIRT}/${image_reference_name}_to_${image_name}_warpcoef_inverted.nii"

#~ # Define labelled reference input paths
#~ image_basalganglia="${folder_atlas_annotated}/AMBMC-c57bl6-basalganglia-labels-15um_reoriented.nii.gz"
#~ image_cerebellum="${folder_atlas_annotated}/AMBMC-c57bl6-cerebellum-labels-15um_reoriented.nii.gz"
#~ image_cortex="${folder_atlas_annotated}/AMBMC-c57bl6-cortex-labels-15um_reoriented.nii.gz"
#~ image_hippocampus="${folder_atlas_annotated}/AMBMC-c57bl6-hippocampus-labels-15um_reoriented.nii.gz"
#~ image_structures=($image_reference $image_basalganglia $image_cerebellum $image_cortex $image_hippocampus)

#~ # Define labelled reference outputs paths which are warped to native space
#~ image_reference_invwarped="${folder_FNIRT}/model_to_${image_name}_invwarped.nii.gz"
#~ image_basalganglia_invwarped="${folder_FNIRT}/basalganglia_to_${image_name}_invwarped.nii"
#~ image_cerebellum_invwarped="${folder_FNIRT}/cerebellum_to_${image_name}_invwarped.nii"
#~ image_cortex_invwarped="${folder_FNIRT}/cortex_to_${image_name}_invwarped.nii"
#~ image_hippocampus_invwarped="${folder_FNIRT}/hippocampus_to_${image_name}_invwarped.nii"
#~ image_structures_invwarped=($image_reference_invwarped $image_basalganglia_invwarped $image_cerebellum_invwarped $image_cortex_invwarped $image_hippocampus_invwarped)

# Define annotation volume paths
image_allen_annotation="${folder_atlas_reference}/annotation_25_to_AMBMC_flirted.nii.gz" # input?
image_allen_annotation_invwarped="${folder_FNIRT}/allen_annotation_to_${image_name}_invwarped.nii" # input?
image_allen_annotation_invflirted="${folder_FNIRT}/allen_annotation_to_${image_name}_invwarped.nii"

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
		-verbose 1



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



#~ # FNIRT reference directly to native space
#~ image_warpaffine_direct="${folder_FNIRT}/${image_name}_to_${image_reference_name}_warpaffine_direct.mat"
#~ image_reference_flirted_direct="${folder_FNIRT}/${image_reference_name}_flirted_direct.nii"
#~ image_warpcoef_direct="${folder_FNIRT}/${image_name}_to_${image_reference_name}_warpcoef_direct.nii"
#~ image_reference_warped_direct="${folder_FNIRT}/${image_reference_name}_to_${image_name}_warped_direct.nii"

#~ echo
#~ echo flirt
#~ flirt 	-in $image_reference \
		#~ -ref $image_inmasked \
		#~ -omat $image_warpaffine_direct \
		#~ -out $image_reference_flirted_direct \
		#~ -verbose 1

#~ echo
#~ echo fnirt
#~ fnirt 	--config=AMBMC_config.cnf \
		#~ --in=$image_reference \
		#~ --ref=$image \
		#~ --cout=$image_warpcoef_direct \
		#~ --aff=$image_warpaffine_direct \
		#~ --refmask=$image_inmask_bin \
		#~ --applyrefmask=1 \
		#~ --iout=$image_reference_warped_direct \
		#~ --verbose
		
#~ # Define labelled reference outputs paths which are warped to native space
#~ image_reference_invwarped_direct="${folder_FNIRT}/model_to_${image_name}_invwarped_direct.nii.gz"
#~ image_basalganglia_invwarped_direct="${folder_FNIRT}/basalganglia_to_${image_name}_invwarped_direct.nii"
#~ image_cerebellum_invwarped_direct="${folder_FNIRT}/cerebellum_to_${image_name}_invwarped_direct.nii"
#~ image_cortex_invwarped_direct="${folder_FNIRT}/cortex_to_${image_name}_invwarped_direct.nii"
#~ image_hippocampus_invwarped_direct="${folder_FNIRT}/hippocampus_to_${image_name}_invwarped_direct.nii"
#~ image_structures_invwarped_direct=($image_reference_invwarped_direct $image_basalganglia_invwarped_direct $image_cerebellum_invwarped_direct $image_cortex_invwarped_direct $image_hippocampus_invwarped_direct)
		
#~ # Warp labelled references to native space
#~ echo
#~ echo "inverted applywarps"
#~ for i in ${!image_structures[@]}; do
	#~ applywarp 	--ref=$image \
				#~ --in=${image_structures[i]} \
				#~ --warp=$image_warpcoef_direct \
				#~ --out=${image_structures_invwarped_direct[i]} \
				#~ --interp=nn \
				#~ --verbose
#~ done
		


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

#~ image_reference_name="${image_reference##*/}"
#~ image_reference_name="${image_reference_name%%.*}"

		

#~ # Warp labelled references to native space
#~ echo
#~ echo "inverted applywarps"
#~ for i in ${!image_structures[@]}; do
	#~ applywarp 	--ref=$image \
				#~ --in=${image_structures[i]} \
				#~ --warp=$image_warpcoef_inverted \
				#~ --out=${image_structures_invwarped[i]} \
				#~ --interp=nn \
				#~ --verbose
#~ done



end=`date +%h`
runtime=$((end-start))
echo $runtime
