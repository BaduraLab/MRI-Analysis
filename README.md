# List of scripts with short explanations
analysis_VOI.py: Compute center of masses for Volumes of Interest and try map volume integers onto a 1D MDS space
BrainExtraction_RATS.sh: Use RATS to generate brain mask images for all files within input folder
Bru2nii_transition.sh: Script to reorient Bru2nii output
cerebellar_adjustment_human.py: Use manually adjusted human cerebellar mask to adjust annotations using nearest neighbour interpolation 
cerebellar_adjustment_mouse.py: Use manually adjusted mouse cerebellar mask to adjust annotations using nearest neighbour interpolation 
compute_mask_volumes.py: Compute volumes of mouse reference images
compute_volumes_human.py: Compute volumes of human subjects (adjusted with manual cerebellar and lobular annotations)
compute_volumes_mouse.py: Compute volumes of human subjects (adjusted with manual cerebellar and lobular annotations)\
folders2gifs.py: Convert folder with 
functions.py: Script that defines helper functions
getNativeCoordinates.py: Compute coordinates in native space from coordinates in reference space for human and mouse subjects
high_detail_to_low_detail.py: Use structure path to change volume integers to their parent volume integers: Use structure path to change volume integers to their parent volume integers
lobular_adjustment_human.py: Use manually adjusted human lobular annotations to adjust cerebellar adjusted annotations using nearest neighbour interpolation
lobular_adjustment_mouse.py: Use manually adjusted mouse lobular annotations to adjust cerebellar adjusted annotations using nearest neighbour interpolation
mask_adjustment_mouse.py: Use RATS generated mask to adjust annotation
process_human.py: Linearly and non-linearly transform human subjects to reference space
process_mouse.py: Linearly and non-linearly transform mouse subjects to reference space