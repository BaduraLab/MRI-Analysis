//run("Script...");
//selectWindow("AVG_control5_reoriented.nii");
//selectWindow("control5_reoriented.nii");
//run("Z Project...");
//run("Z Project...", "projection=[Max Intensity]");
//selectWindow("control5_reoriented.nii");
//selectWindow("AVG_control5_reoriented.nii");
//selectWindow("MAX_control5_reoriented.nii");


// loop over structures, control4 and patient, adjust intensity scale to max(control4, patient) use annotation files only and use sum Zproj (generally best way to show IMO)


controlIsolationsFolder = "/mnt/tosh/Projects/MEP/mep-scripts/Data/Human/Processed/control4/isolations"
patientIsolationsFolder = "/mnt/tosh/Projects/MEP/mep-scripts/Data/Human/Processed/patient/isolations"

print(controlIsolationsFolder);
print(patientIsolationsFolder);

//controlIsolationsFiles = getFileList(controlIsolationsFolder + File.separator + "*" + File.separator + "*annotation*_COMcentered.nii.gz");
//patientIsolationsFiles = getFileList(patientIsolationsFolder + File.separator + "*" + File.separator + "*annotation*_COMcentered.nii.gz");
controlIsolationsFile = controlIsolationsFolder + File.separator + "si" + File.separator + "control4_annotation_subcortical_thrarg_flirtedRigid_ventral tegmental area_COMcentered.nii.gz";
patientIsolationsFile = patientIsolationsFolder + File.separator + "si" + File.separator + "patient_annotation_subcortical_thrarg_flirtedRigid_ventral tegmental area_COMcentered.nii.gz";

print(controlIsolationsFile);
print(patientIsolationsFile);

//for(i=0;i<controlIsolationsFiles.length;i++){
//	controlIsolationsFile = controlIsolationsFiles[i];
//	patientIsolationsFile = patientIsolationsFiles[i];
//	print(controlIsolationsFile);
//	print(patientIsolationsFile);
//}

print("open=[" + controlIsolationsFile + "]")

run("NIfTI-Analyze", "open=[" + controlIsolationsFile + "]");
run("NIfTI-Analyze", "open=" + patientIsolationsFile);

//getFileList(controlIsolationsFolder);

//files = getFileList(dir+dirs[u]+"processed"+File.separator+"transformations"+File.separator);
//print(files.length);
//for(i=0;i<2.1;i++){
//for(i=0;i<files.length;i++){

//dir = getDirectory("Choose Data Directory");
//print(dir);
//main = dir+dirs[u]+"processed"+File.separator+"transformations"+File.separator+files[i];
//open(main);
//run("Z Project...", "projection=[Sum Intensity]");
//save(path);
//run("Close All");

// get files in isolation folder
// run Zprojection
// save folder
