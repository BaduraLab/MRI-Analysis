raw_folder="/home/enzo/Desktop/Data/Mouse/Raw"

custom_array=("${raw_folder}/KO_60_female.nii")
for subject in "${custom_array[@]}" ; do

#~ subject="WT_50.nii"
echo ${subject}
subject_name=${subject##*/}
subject_name=${subject_name%.*}
echo $subject_name
raw="${subject}"

processed_new_folder="/home/enzo/Desktop/Data/Mouse/Processed_New"
mkdir -p $processed_new_folder
processed_folder="/home/enzo/Desktop/Data/Mouse/Processed/${subject_name}"
processed_new_subject_folder="${processed_new_folder}/${subject_name}"
mkdir -p $processed_new_subject_folder
processed_new="${processed_new_subject_folder}/${subject_name}.nii.gz"

processed_mask="${processed_folder}/${subject_name}_mask_t=500_v=380_k=6.nii"
processed_vtp="${processed_folder}/${subject_name}_vtp_t=500_v=380_k=6_a=0_n=2000.vtp"
processed_vtp_mask="${processed_folder}/${subject_name}_vtp_mask_t=500_v=380_k=6_a=0_n=2000.nii"
processed_mask_manual="${processed_folder}/${subject_name}_mask_t=500_v=380_k=6.mask.nii"

processed_new_mask="${processed_new_subject_folder}/${subject_name}_mask_t=500_v=380_k=6.nii.gz"
processed_new_vtp="${processed_new_subject_folder}/${subject_name}_vtp_t=500_v=380_k=6_a=0_n=2000.vtp"
processed_new_vtp_mask="${processed_new_subject_folder}/${subject_name}_vtp_mask_t=500_v=380_k=6_a=0_n=2000.nii.gz"
processed_new_mask_manual="${processed_new_subject_folder}/${subject_name}_mask_t=500_v=380_k=6.mask.nii.gz"

gzip -c $raw > $processed_new
gzip -c $processed_mask > $processed_new_mask
cp $processed_vtp $processed_new_vtp
gzip -c $processed_vtp_mask > $processed_new_vtp_mask
gzip -c $processed_mask_manual > $processed_new_mask_manual

fslorient -copysform2qform $processed_new
fslorient -setqformcode 1 $processed_new
fslreorient2std $processed_new $processed_new

masks=($processed_new_mask $processed_new_mask_vtp $processed_new_mask_manual)
for mask in "${masks[@]}"; do
	fslcpgeom $processed_new $mask
	fslswapdim $mask x y z $mask
done

fsleyes $processed_new $processed_new_mask $processed_new_vtp_mask $processed_new_mask_manual

done
