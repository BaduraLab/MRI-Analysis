# Define paths
atlas_folder=/usr/local/fsl/data/standard
allen_annotation=${atlas_folder}/allen/annotation_25_reoriented.nii.gz
allen_average_template=${atlas_folder}/allen/average_template_25_reoriented.nii.gz
AMBMC_average_template=${atlas_folder}/AMBMC_model_reoriented.nii.gz

allen_to_AMBMC_flirt=${atlas_folder}/allen/average_template_25_to_AMBMC_flirt.mat
allen_average_template_to_AMBMC_flirted=${atlas_folder}/allen/average_template_25_to_AMBMC_flirt.nii.gz
allen_annotation_to_AMBMC_flirted=${atlas_folder}/allen/annotation_25_to_AMBMC_flirted.nii.gz

allen_to_AMBMC_fnirt=${atlas_folder}/allen/average_template_25_to_AMBMC_fnirt.nii.gz
allen_average_template_to_AMBMC_fnirted=${atlas_folder}/allen/average_template_25_to_AMBMC_fnirted.nii.gz
allen_annotation_to_AMBMC_fnirted=${atlas_folder}/allen/annotation_25_to_AMBMC_fnirted.nii.gz

allen_annotation_bin=${atlas_folder}/allen/annotation_25_reoriented_bin.nii.gz
allen_annotation_bin_flirted=${atlas_folder}/allen/annotation_25_reoriented_bin_flirted.nii.gz



# Create Allen mask
fslmaths $allen_annotation -thr 0.5 -bin $allen_annotation_bin



# FLIRT (native2reference)
echo
echo flirt
flirt 	-in $allen_average_template \
		-ref $AMBMC_average_template \
		-omat $allen_to_AMBMC_flirt \
		-out $allen_average_template_to_AMBMC_flirted \
		-verbose 1

# FLIRT applywarp annotation (native2reference)
echo
echo "Invert FLIRT allen to image"
flirt 	-in $allen_annotation \
		-ref $AMBMC_average_template \
		-out $allen_annotation_to_AMBMC_flirted \
		-init $allen_to_AMBMC_flirt \
		-applyxfm \
		-verbose 1
		
# FLIRT applywarp annotation (native2reference)
echo
echo "Invert FLIRT allen to image"
flirt 	-in $allen_annotation_bin \
		-ref $AMBMC_average_template \
		-out $allen_annotation_bin_to_AMBMC_flirted \
		-init $allen_to_AMBMC_flirt \
		-applyxfm \
		-verbose 1
		
		
		
# FNIRT (native2reference)
echo
echo fnirt
fnirt 	--config=AMBMC_config.cnf \
		--in=$allen_average_template \
		--ref=$AMBMC_average_template \
		--cout=$allen_to_AMBMC_fnirt \
		--aff=$allen_to_AMBMC_flirt \
		--iout=$allen_average_template_to_AMBMC_fnirted \
		--verbose
		
# FNIRT applywarp annotation (native2reference)
echo
echo applywarp
applywarp 	--ref=$AMBMC_average_template \
			--in=$allen_annotation \
			--warp=$allen_to_AMBMC_fnirt \
			--out=$allen_annotation_to_AMBMC_fnirted \
			--interp=nn \
			--verbose
