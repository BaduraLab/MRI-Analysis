# Decompose mouse flirt affine and wrtie flirt rigid as transformation component
data_folder = file.path(getwd(), "Data", "Mouse", "Processed_Old")
Sys.glob(file.path(data_folder, "*"), dirmark=TRUE)
mouse_folders = list.dirs(data_folder, recursive = FALSE)
for (mouse_folder in mouse_folders) {
  print(mouse_folder)
  affine_path = file.path(mouse_folder, "flirtAffine.mat")
  flirtAffine = read.delim(affine_path, header=FALSE, sep = " ")
  flirtAffine = RNiftyReg::asAffine(data.matrix(flirtAffine[,c(1,3,5,7)]))
  flirtAffine_decomposed = RNiftyReg::decomposeAffine(flirtAffine)
  flirtRigid = diag(4)
  flirtRigid[1:3, 1:3] = flirtAffine_decomposed$rotationMatrix
  #flirtRigid[1:3, 4] = flirtAffine_decomposed$translation
  
  write.table(flirtRigid, file=file.path(mouse_folder, "flirtRigid.mat"), row.names=FALSE, col.names=FALSE, sep = "  ")
}



# loop through files
# decompose affine
# compose rigid transformation matrx out of rotation and translation parameters
# save rigid transformation matrix



# Decompose human flirt affine and wrtie flirt rigid as transformation component
data_folder = file.path(getwd(), "Data", "Human", "Processed")
Sys.glob(file.path(data_folder, "*"), dirmark=TRUE)
subject_folders = list.dirs(data_folder, recursive = FALSE)
for (subject_folder in subject_folders) {
  print(subject_folder)
  affine_path = file.path(subject_folder, "flirtAffine.mat")
  flirtAffine = read.delim(affine_path, header=FALSE, sep = " ")
  flirtAffine = RNiftyReg::asAffine(data.matrix(flirtAffine[,c(1,3,5,7)]))
  flirtAffine_decomposed = RNiftyReg::decomposeAffine(flirtAffine)
  flirtRigid = diag(4)
  flirtRigid[1:3, 1:3] = flirtAffine_decomposed$rotationMatrix
  #flirtRigid[1:3, 4] = flirtAffine_decomposed$translation
  flirtRigid[1:3, 4] = flirtAffine[1:3, 4]
  
  write.table(flirtRigid, file=file.path(subject_folder, "flirtRigid.mat"), row.names=FALSE, col.names=FALSE, sep = "  ")
}
