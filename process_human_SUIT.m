%suit
t1_path = {'C:\Users\Enzon\Documents\Projects\MEP\mep-scripts\Data\Human\Processed\patient\patient_reoriented.nii', 'C:\Users\Enzon\Documents\Projects\MEP\mep-scripts\Data\Human\Processed\control2\control2_reoriented.nii', 'C:\Users\Enzon\Documents\Projects\MEP\mep-scripts\Data\Human\Processed\control1\control1_reoriented.nii'};
% t1_path{1}
suit_isolate_seg(t1_path(1))
suit_isolate_seg(t1_path(2))
suit_isolate_seg(t1_path(3))

nT1 = length(t1_path);
for iT1 = 1:nT1
    [t1_filepath, t1_filename, ~] = fileparts(t1_path{iT1});
    
    job.subjND.gray = {[t1_filepath, filesep, t1_filename, '_seg1.nii']};
    job.subjND.white = {[t1_filepath, filesep, t1_filename, '_seg2.nii']};
    job.subjND.isolation = {[t1_filepath, filesep, 'c_', t1_filename, '_pcereb.nii']};
    
    suit_normalize_dartel(job)
end

clear job
V = cell(1, nT1);
for iT1 = 1:nT1
    [t1_filepath, t1_filename, ~] = fileparts(t1_path{iT1});
    
    job.Affine = {[t1_filepath, filesep, 'Affine_', t1_filename, '_seg1.mat']};
    job.flowfield = {[t1_filepath, filesep, 'u_a_', t1_filename, '_seg1.nii']};
    job.ref = t1_path(iT1);
    job.resample = {'C:\Program Files\MATLAB\R2020a\toolbox\spm12\toolbox\suit\atlasesSUIT\Lobules-SUIT.nii'};
    
    suit_reslice_dartel_inv(job)
    
    V{iT1} = suit_vol([t1_filepath, filesep, 'iw_Lobules-SUIT_u_a_', t1_filename, '_seg1.nii']);
end