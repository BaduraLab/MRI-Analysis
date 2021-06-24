#!/bin/bash

# Define package folder
package_folder="/home/enzo/Downloads/ambmc-c57bl6-FSL-atlas_v0.8"

# Move files to appropriate folders
target_folder="${FSLDIR}/data/standard/"
sudo cp -p $package_folder/AMBMC_model.nii.gz $target_folder

sudo mkdir $FSLDIR/data/atlases/AMBMC/
target_folder="${FSLDIR}/data/atlases/AMBMC/" 
sudo cp -p "${package_folder}/AMBMC_model.nii.gz" $target_folder
sudo cp -p $package_folder/AMBMC-c57bl6-cerebellum-labels-15um.nii.gz $target_folder
sudo cp -p $package_folder/AMBMC-c57bl6-cortex-labels-15um.nii.gz $target_folder
sudo cp -p $package_folder/AMBMC-c57bl6-hippocampus-labels-15um.nii.gz $target_folder

target_folder="$FSLDIR/data/atlases/" 
sudo cp -p $package_folder/AMBMC_basalganglia_labels.xml $target_folder
sudo cp -p $package_folder/AMBMC_cerebellum_labels.xml $target_folder
sudo cp -p $package_folder/AMBMC_cortex_labels.xml $target_folder
sudo cp -p $package_folder/AMBMC_hippocampus_labels.xml $target_folder

target_folder="$FSLDIR/etc/luts/"
sudo cp -p $package_folder/ambmc_cortex_LUT.rgb $target_folder
sudo cp -p $package_folder/ambmc_basal_LUT.rgb $target_folder
sudo cp -p $package_folder/ambmc_cereb_LUT.rgb $target_folder

target_folder="$FSLDIR/etc/flirtsch/"
sudo cp -p $package_folder/AMBMC_config.cnf $target_folder
