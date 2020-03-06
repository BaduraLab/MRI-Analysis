#!/bin/bash

# Define folder and file paths
folder_working="/home/enzo/Desktop" # Alternatively define this by current working directory
folder_atlas_reference="${FSLDIR}/data/standard"
folder_atlas_annotated="${FSLDIR}/data/atlases/AMBMC/"
folder_atlas_xml="${FSLDIR}/data/atlases/"
folder_atlas_LUT="${FSLDIR}/etc/luts/"
folder_atlas_config="${FSLDIR}/etc/flirtsch/"

image="${folder_working}/Data/Mouse/Test/KO_2.nii" #~ input?
image_dir="${image%.*}"
mkdir -p $image_dir
image_name="${image##*/}"

image_reference="${folder_atlas_reference}/AMBMC_model.nii.gz" #~ input?
image_reference_name="${image_reference##*/}"

image_warpaff="${image_dir}/${image_name%.*}_to_${image_reference_name%.*}_warpaff.mat"
image_warpcoef="${image_dir}/${image_name%.*}_to_${image_reference_name%.*}_warpcoef.nii"
image_warped="${image_dir}/${image_name%.*}_to_${image_reference_name%.*}_warped.nii"
image_warpcoef_inverted="${image_dir}/${image_reference_name%.*}_to_${image_name%.*}_warpcoef_inverted.nii"

image_basalganglia="${folder_atlas_annotated}/AMBMC-c57bl6-basalganglia-labels-15um.nii.gz"
image_cerebellum="${folder_atlas_annotated}/AMBMC-c57bl6-cerebellum-labels-15um.nii.gz"
image_cortex="${folder_atlas_annotated}/AMBMC-c57bl6-cortex-labels-15um.nii.gz"
image_hippocampus="${folder_atlas_annotated}/AMBMC-c57bl6-hippocampus-labels-15um.nii.gz"
image_structures=(image_basalganglia image_cerebellum image_cortex image_hippocampus)

image_basalganglia_invwarped="${image_dir}/basalganglia_to_${image_name%.*}_invwarped.nii"
image_cerebellum_invwarped="${image_dir}/cerebellum_to_${image_name%.*}_invwarped.nii"
image_cortex_invwarped="${image_dir}/cortex_to_${image_name%.*}_invwarped.nii"
image_hippocampus_invwarped="${image_dir}/hippocampus_to_${image_name%.*}_invwarped.nii"
image_structures_invwarped=(image_basalganglia_invwarped image_cerebellum_invwarped image_cortex_invwarped image_hippocampus_invwarped)

#~ filename=$(basename -- "$fullfile")
#~ extension="${filename##*.}"
#~ filename="${filename%.*}"
#~ filename="${fullfile##*/}"

#~ # Move files to appropriate folders
#~ target_folder="${FSLDIR}/data/standard/"
#~ sudo cp -p $package_folder/AMBMC_model.nii.gz $target_folder

#~ sudo mkdir ${FSLDIR}/data/atlases/AMBMC/
#~ target_folder="${FSLDIR}/data/atlases/AMBMC/" 
#~ sudo cp -p "${package_folder}/AMBMC-c57bl6-basalganglia-labels-15um.nii.gz" $target_folder
#~ sudo cp -p $package_folder/AMBMC-c57bl6-cerebellum-labels-15um.nii.gz $target_folder
#~ sudo cp -p $package_folder/AMBMC-c57bl6-cortex-labels-15um.nii.gz $target_folder
#~ sudo cp -p $package_folder/AMBMC-c57bl6-hippocampus-labels-15um.nii.gz $target_folder

#~ target_folder="${FSLDIR}/data/atlases/" 
#~ sudo cp -p $package_folder/AMBMC_basalganglia_labels.xml $target_folder
#~ sudo cp -p $package_folder/AMBMC_cerebellum_labels.xml $target_folder
#~ sudo cp -p $package_folder/AMBMC_cortex_labels.xml $target_folder
#~ sudo cp -p $package_folder/AMBMC_hippocampus_labels.xml $target_folder

#~ target_folder="${FSLDIR}/etc/luts/"
#~ sudo cp -p $package_folder/ambmc_cortex_LUT.rgb $target_folder
#~ sudo cp -p $package_folder/ambmc_basal_LUT.rgb $target_folder
#~ sudo cp -p $package_folder/ambmc_cereb_LUT.rgb $target_folder

#~ target_folder="${FSLDIR}/etc/flirtsch/"
#~ sudo cp -p $package_folder/AMBMC_config.cnf $target_folder

flirt -ref $image_reference -in $image -omat $image_warpaff

applywarp --ref=$image_reference --in=$image --out=$image_warped --warp=$image_warpaff

#~ fnirt 	--config=AMBMC_config.cnf \
		#~ --in=$image \
		#~ --ref=$image_reference \
		#~ --cout=$image_warpcoef \
		#~ --iout=$image_warped \
		#~ --aff=$image_warp

#~ applywarp 	--in=/home/enzo/Desktop/wetransfer-489aa5/Pax5_WT/WT_10b.nii \
			#~ --ref=/usr/local/fsl/data/standard/AMBMC_model.nii.gz \
			#~ --warp=/home/enzo/Desktop/wetransfer-489aa5/Pax5_WT/WT_10b_warpcoef.nii.gz \
			#~ --out=/home/enzo/Desktop/wetransfer-489aa5/Pax5_WT/WT_10b_warped.nii \
			/

#~ invwarp --ref=$image \
		#~ --warp=$image_warpcoef \
		#~ --out=$image_warpcoef_inverted
		
#~ for i in ${!image_structures[@]}; do
	#~ applywarp 	--ref=$image \
				#~ --in=${image_structures[i]} \
				#~ --warp=image_warpcoef_inverted \
				#~ --out=${image_structures_invwarped[i]}
#~ done
