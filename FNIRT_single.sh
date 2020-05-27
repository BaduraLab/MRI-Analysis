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

subject_name="WT_50"
folder_working="/home/enzo/Desktop/Data/Mouse/Processed_New/${subject_name}" # Alternatively define this by current working directory
folder_FLIRT="${folder_working}/FLIRT"
datestr="$(date +%F)-$(date +%H)-$(date +%M)-$(date +%S)"
datestr="2020-05-17-10-49-04" ##############################################

folder_FNIRT="${folder_working}/FNIRT_${datestr}"
mkdir -p $folder_FNIRT
cp /home/enzo/Desktop/mep-scripts/FNIRT_single.sh $folder_FNIRT/

# Log files
FNIRT_initial_log="${folder_FNIRT}/FNIRT_initial_log"
FNIRT_log="${folder_FNIRT}/FNIRT_log"

# Folders of reference
folder_atlas_reference="${FSLDIR}/data/standard/allen_new"
folder_atlas_annotated="${FSLDIR}/data/atlases/AMBMC"
folder_atlas_xml="${FSLDIR}/data/atlases"
folder_atlas_LUT="${FSLDIR}/etc/luts"
folder_atlas_config="${FSLDIR}/etc/flirtsch"

# Define image path and reorient, use reoriented image as image
image="${folder_FLIRT}/${subject_name}.nii" # input?

# Define image mask path and reorient, use reoriented mask as mask
image_inmask="${folder_FLIRT}/${subject_name}_mask_t=500_v=380_k=6.mask.nii.gz"

# Define image and mask FLIRT output paths
image_inmask_bin="${folder_FLIRT}/${subject_name}_mask_t=500_v=380_k=6_bin.mask.nii.gz"
image_inmask_bin_warped="${folder_FNIRT}/${subject_name}_mask_t=500_v=380_k=6_bin_warped.mask.nii.gz"
image_inmask_bin_warpcoef="${folder_FNIRT}/${subject_name}_mask_t=500_v=380_k=6_bin_warpcoef.mask.nii.gz"

image_name="${image##*/}"
image_name="${image_name%.*}"

# Define image reference path and check whether it is reoriented to reference
image_reference="${folder_atlas_reference}/average_template_25_to_AMBMC_flirted.nii.gz" # input?
image_reference_mask="${folder_atlas_reference}/annotation_25_bin_to_AMBMC_flirted.nii.gz" # input?
image_reference_name="allen_model" # input? average_template_25_to_AMBMC_flirted
echo "Print reorientation-of-reference-to-standard-reference matrix"
fslreorient2std $image_reference # Prints matrix, if identity then ok

# Define image warped with warp between in and ref masks
image_warped_torefmask="${folder_FNIRT}/${image_name}_to_${image_reference_name}_warped_torefmask.nii"
image_allen_annotation_torefmask_invwarped="${folder_FNIRT}/allen_annotation_warped_torefmask.nii"

# Define FLIRT and FNIRT paths
image_warpaffine="${folder_FLIRT}/${image_name}_to_${image_reference_name}_warpaffine.mat"
image_warpcoef="${folder_FNIRT}/${image_name}_to_${image_reference_name}_warpcoef.nii"
image_warped="${folder_FNIRT}/${image_name}_to_${image_reference_name}_warped.nii"
image_warpaffine_inverted="${folder_FLIRT}/${image_reference_name}_to_${image_name}_warpaffine_inverted.mat"
image_warpcoef_inverted="${folder_FNIRT}/${image_reference_name}_to_${image_name}_warpcoef_inverted.nii"
image_inmask_bin_warpcoef_inverted="${folder_FNIRT}/${subject_name}_mask_t=500_v=380_k=6_bin_warpcoef_inverted.mask.nii.gz"

# Define annotation volume paths
image_allen_annotation="${folder_atlas_reference}/annotation_25_to_AMBMC_flirted.nii.gz" # input?
image_allen_annotation_bin="${folder_atlas_reference}/annotation_25_reoriented_bin_to_AMBMC_flirted.nii.gz"
image_allen_annotation_invwarped="${folder_FNIRT}/allen_annotation_to_${image_name}_invwarped.nii" # input?
image_allen_annotation_invflirted="${folder_FLIRT}/allen_annotation_to_${image_name}_invflirted.nii"



## FNIRT initial guess by warping mask to reference(native2reference)
#echo "fnirt initial"
#both $FNIRT_initial_log \
#fnirt 	--in=$image_inmask_bin \
		#--ref=$image_reference_mask \
		#--cout=$image_inmask_bin_warpcoef \
		#--iout=$image_inmask_bin_warped \
		#--applyinmask=0 \
		#--applyrefmask=0 \
		#--imprefm=0 \
		#--impinm=0 \
		#--verbose \
		#--lambda=150,50,20,10,5,1,0.1 \
		#--infwhm=0.18,0.15,0.12,0.09,0.06,0,0 \
		#--reffwhm=0.06,0.06,0.03,0.03,0,0,0 \
		#--miter=5,5,5,5,5,10,10 \
		#--subsamp=8,8,4,4,4,4,2 \
		#--estint=0 \
		#--ssqlambda=1 \
		#--warpres=0.66,1,0.46 \
		#--regmod=bending_energy \
		#--intmod=global_non_linear \
		#--intorder=5 \
		#--biasres=3,5,2 \
		#--biaslambda=10000 \
		#--refderiv=0 \
		#--numprec=float \
		#--splineorder=3 \
		#--aff=$image_warpaffine

## FNIRT applywarp (native2refmask)
#echo
#echo "apply initial warp to the reference mask to image"
#applywarp 	--ref=$image_reference_mask \
			#--in=$image \
			#--warp=$image_inmask_bin_warpcoef \
			#--out=$image_warped_torefmask \
			#--verbose



## FNIRT (native2reference) parameter discrepancy test
## Currently parameter differences are reduced to applyinmask=0->1 and estint=0->1
## Of course in, ref, cout and iout are changed appropriately
## toorig: subsamp, miter, in/reffwhm, lambda, estint
#echo "fnirt discrepancy test"
#both $FNIRT_log fnirt \
	#--in=$image_warped_torefmask \
	#--ref=$image_reference \
	#--cout=$image_warpcoef \
	#--iout=$image_warped \
	#--inmask=$image_inmask_bin_warped \
	#--applyinmask=1 \
	#--applyrefmask=0 \
	#--imprefm=0 --impinm=0 \
	#--verbose \
	#--subsamp=8 \
	#--miter=15 \
	#--infwhm=0.18 \
	#--reffwhm=0.06 \
	#--lambda=150 \
	#--estint=1 \
	#--ssqlambda=0 \
	#--warpres=0.66,1,0.46 \
	#--regmod=bending_energy \
	#--intmod=global_non_linear \
	#--intorder=5 \
	#--biasres=3,5,2 \
	#--biaslambda=10000 \
	#--refderiv=0 \
	#--numprec=float \
	#--splineorder=3
##~ fnirt --in=/home/enzo/Desktop/Data/Mouse/Processed_New/WT_50/FNIRT_2020-05-06-13-02-13/WT_50_to_allen_model_warped_torefmask.nii --ref=/usr/local/fsl/data/standard/allen/average_template_25_to_AMBMC_flirted.nii.gz --cout=/home/enzo/Desktop/Data/Mouse/Processed_New/WT_50/FNIRT_2020-05-06-13-02-13/WT_50_to_allen_model_warpcoef.nii --iout=/home/enzo/Desktop/Data/Mouse/Processed_New/WT_50/FNIRT_2020-05-06-13-02-13/WT_50_to_allen_model_warped.nii --inmask=/home/enzo/Desktop/Data/Mouse/Processed_New/WT_50/FNIRT_2020-05-06-13-02-13/WT_50_mask_t=500_v=380_k=6_bin_warped.mask.nii.gz --applyinmask=1 --applyrefmask=0 --imprefm=0 --impinm=0 --verbose --subsamp=4 --miter=30 --infwhm=0.18 --reffwhm=0.06 --lambda=150 --estint=1 --ssqlambda=0 --warpres=0.66,1,0.46 --regmod=bending_energy --intmod=global_non_linear --intorder=5 --biasres=3,5,2 --biaslambda=10000 --refderiv=0 --numprec=float --splineorder=3	
##~ echo "fnirt discrepancy test"
##~ both $FNIRT_log fnirt --in=$image --ref=$image_reference --cout=$image_warpcoef --iout=$image_warped --inwarp=$image_inmask_bin_warpcoef --inmask=$image_inmask_bin --applyinmask=0 --applyrefmask=0 --imprefm=0 --impinm=0 --verbose --subsamp=4 --miter=30 --infwhm=0.18 --reffwhm=0.06 --lambda=150 --estint=1 --ssqlambda=0 --warpres=0.66,1,0.46 --regmod=bending_energy --intmod=global_non_linear --intorder=5 --biasres=3,5,2 --biaslambda=10000 --refderiv=0 --numprec=float --splineorder=3



## Invert FNIRT (reference2native)
#echo
#echo "Invert FNIRT"
#invwarp --ref=$image \
		#--warp=$image_inmask_bin_warpcoef \
		#--out=$image_inmask_bin_warpcoef_inverted \
		#--verbose
		
# Invert FNIRT (reference2native)
echo
echo "Invert FNIRT"
invwarp --ref=$image_warped_torefmask \
		--warp=$image_warpcoef \
		--out=$image_warpcoef_inverted \
		--verbose



# FNIRT applywarp (reference2native)
# Warp labelled references to native space
echo
echo "inverted applywarp"
applywarp 	--ref=$image_warped_torefmask \
			--in=$image_allen_annotation \
			--warp=$image_warpcoef_inverted \
			--out=$image_allen_annotation_torefmask_invwarped \
			--interp=nn \
			--verbose

# FNIRT applywarp (reference2native)
# Warp labelled references to native space
echo
echo "inverted applywarp"
applywarp 	--ref=$image \
			--in=$image_allen_annotation_torefmask_invwarped \
			--warp=$image_inmask_bin_warpcoef_inverted \
			--out=$image_allen_annotation_invwarped \
			--interp=nn \
			--verbose



end=`date +%h`
runtime=$((end-start))
echo $runtime



		
#~ #	name of reference image
#~ --ref=/usr/local/fsl/data/standard/AMBMC_model.nii.gz
#~ #	If =1, use implicit masking based on value in --ref image. Default =1
#~ --imprefm=1
#~ #	If =1, use implicit masking based on value in --in image, Default =1
#~ --impinm=1
#~ --impinval=0
#~ #	sub-sampling scheme
#~ --subsamp=4,4,2,2,1,1
#~ # 	Max # of non-linear iterations
#~ --miter=5,5,5,5,5,10
#~ #	FWHM (in mm) of gaussian smoothing kernel for input volume
#~ --infwhm=0.18,0.15,0.12,0.09,0.06,0
#~ #	FWHM (in mm) of gaussian smoothing kernel for ref volume
#~ --reffwhm=0.06,0.06,0.03,0.03,0,0
#~ #	Weigth of membrane energy regularisation, default depending on --ssqlambda and --regmod switches. See user documetation.
#~ --lambda=150,75,50,30,20,15
#~ #	Estimate intensity-mapping if set, deafult 1 (true)
#~ --estint=1,1,1,1,1,0
#~ #       Apply the mask if set, default 1 (true)
#~ --applyrefmask=0,0,0,0,0,0
#~ #       Apply the mask if set, default 1 (true)
#~ --applyinmask=0
#~ #	(approximate) resolution (in mm) of warp basis in x-, y- and z-direction
#~ --warpres=0.66,1,0.46
#~ #	If set (=1), lambda is weighted by current ssq, default 1
#~ --ssqlambda=1
#~ #	Model for regularisation of warp-field [membrane_energy bending_energy], default bending_energy
#~ --regmod=bending_energy
#~ #	Model for intensity-mapping [none global_linear global_non_linear local_linear global_non_linear_with_bias local_non_linear]
#~ --intmod=global_non_linear
#~ #	Order of poynomial for mapping intensities, default 5
#~ --intorder=5
#~ #	Resolution (in mm) of bias-field modelling local intensities. Determines the knot-spacing for the splines that are used to model a bias-field. It means the same thing as --warpres, but for the bias-field rather than the warp-fields. It is relevant for the local_linear, the global_non_linear_with_bias and the local_non_linear models. Typically a bias-field varies quite slowly over space.
#~ --biasres=3,5,2
#~ #	Weight of regularisation for bias-field, default 10000
#~ --biaslambda=10000
#~ #	If =1, ref image is used to calculate derivatives. Default =0
#~ --refderiv=0
#~ # Its value can be either float or double (default) and it specifies the precision that the hessian H is calculated and stored in. Changing this to float will decrease the amount of RAM needed to store H and will hence allow one to go to slightly higher warp-resolution.
#~ --numprec=float
#~ #	Specifies the order of the B-spline functions modelling the warp-fields, default 3
#~ --splineorder=3

#~ both $FNIRT_log \
#~ fnirt 	--in=$image \
		#~ --ref=$image_reference \
		#~ --cout=$image_warpcoef \
		#~ --iout=$image_warped \
		#~ --inwarp=$image_inmask_bin_warpcoef \
		#~ --inmask=$image_inmask_bin \
		#~ --applyinmask=0 \
		#~ --applyrefmask=0 \
		#~ --imprefm=0 \
		#~ --impinm=0 \
		#~ --verbose \
		#~ --subsamp=4 \
		#~ --miter=30 \
		#~ --infwhm=0.18 \
		#~ --reffwhm=0.06 \
		#~ --lambda=150 \
		#~ --estint=1 \
		#~ --ssqlambda=0 \
		#~ --warpres=0.66,1,0.46 \
		#~ --regmod=bending_energy \
		#~ --intmod=global_non_linear \
		#~ --intorder=5 \
		#~ --biasres=3,5,2 \
		#~ --biaslambda=10000 \
		#~ --refderiv=0 \
		#~ --numprec=float \
		#~ --splineorder=3
		
		#~ # FNIRT (native2reference)
#~ echo
#~ echo fnirt
#~ both $FNIRT_log \
#~ fnirt 	--config=AMBMC_config.cnf \
		#~ --in=$image \
		#~ --ref=$image_reference \
		#~ --cout=$image_warpcoef \
		#~ --inwarp=$image_inmask_bin_warpcoef \
		#~ --inmask=$image_inmask_bin \
		#~ --applyinmask=1,1,1,1,1,1 \
		#~ --iout=$image_warped \
		#~ --verbose \
		#~ --subsamp=4,4,2,2,1,1 \
		#~ --miter=5,5,5,5,5,10 \
		#~ --infwhm=0.18,0.15,0.12,0.09,0.06,0 \
		#~ --reffwhm=0.06,0.06,0.03,0.03,0,0 \
		#~ --lambda=150,50,30,15,5,1 \
		#~ --estint=1,1,1,1,1,0 \
		#~ --applyrefmask=0,0,0,0,0,0 \
		#~ --ssqlambda=0
	
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
