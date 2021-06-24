run("NIfTI-Analyze", "open=[/mnt/tosh/Projects/MEP/mep-scripts/Data/Human/Processed/control4/Isolations/si/control4_annotation_subcortical_thrarg_flirtedRigid_ventral tegmental area_COMcentered.nii.gz]");
run("NIfTI-Analyze", "open=[/mnt/tosh/Projects/MEP/mep-scripts/Data/Human/Processed/patient/Isolations/si/patient_annotation_subcortical_thrarg_flirtedRigid_ventral tegmental area_COMcentered.nii.gz]");

selectWindow("control4_annotation_subcortical_thrarg_flirtedRigid_ventral tegmental area_COMcentered.nii.gz");
run("Z Project...", "projection=[Sum Slices]");
selectWindow("SUM_control4_annotation_subcortical_thrarg_flirtedRigid_ventral tegmental area_COMcentered.nii.gz");
save("/mnt/tosh/Projects/MEP/mep-scripts/Data/Human/Processed/control4/Isolations/si/control4_annotation_subcortical_thrarg_flirtedRigid_ventral tegmental area_COMcentered.png");

selectWindow("patient_annotation_subcortical_thrarg_flirtedRigid_ventral tegmental area_COMcentered.nii.gz");
run("Z Project...", "projection=[Sum Slices]");
selectWindow("SUM_patient_annotation_subcortical_thrarg_flirtedRigid_ventral tegmental area_COMcentered.nii.gz");
run("Brightness/Contrast...");
setMinAndMax(0, 66);
save("/mnt/tosh/Projects/MEP/mep-scripts/Data/Human/Processed/patient/Isolations/si/patient_annotation_subcortical_thrarg_flirtedRigid_ventral tegmental area_COMcentered.png");

run("Close All");