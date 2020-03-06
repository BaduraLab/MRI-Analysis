#!/bin/bash

imagename="/home/enzo/Desktop/Data/Mouse/Test/WT_50_vtp_mask_t=500_v=380_k=6_a=0_n=2000.nii"
dir_out="${imagename%.*}.reorient"
image="${dir_out}/${imagename##*/}"

mkdir -p $dir_out

#~ echo $imagename
#~ echo $image
#~ which cp
cp $imagename $image

x="x"
y="y"
z="z"
qformcode="1"

image_deleteorient="${image%.*}_deleterorient.nii"
image_swapdim="${image%.*}_${x}_${y}_${z}.nii.gz"
image_orient="${image%.*}_qformcode=${qformcode}.nii"

echo "Start reorientation"

#~ echo $image
#~ echo $image_deleteorient
#~ cp $image $image_deleteorient
#~ fslorient -deleteorient $image_deleteorient 
#~ fslswapdim $image_deleteorient ${x} ${y} ${z} $image_swapdim
#~ echo $image_swapdim
#~ echo $image_orient
#~ cp $image_swapdim $image_orient
#~ fslorient -setqformcode ${qformcode} $image_orient

#~ fslswapdim $image_orient LR AP IS $image_orient
#~ fslorient -setqformcode ${qformcode} $image_orient

cp $image $image_deleteorient
fslorient -deleteorient $image_deleteorient 
#~ fslswapdim $image_deleteorient ${x} ${y} ${z} $image_swapdim
#~ cp $image_swapdim $image_orient
cp $image_deleteorient $image_orient
fslorient -setqformcode ${qformcode} $image_orient
fslorient -setqform 0.05 0 0 0 0 -0.05 0 0 0 0 0.05 0 0 0 0 1 $image_orient
#~ fslorient -copyqform2sform $image_orient
#~ fslorient -setsform 1 0 0 0 0 -1 0 0 0 0 1 0 0 0 0 1 $image_orient

fslhd $image_orient
fsleyes $image_orient
