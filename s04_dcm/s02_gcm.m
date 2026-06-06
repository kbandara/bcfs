%Collate DCMs into a GCM 
%kav 2025

addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12');
spm('defaults', 'eeg');
spm_jobman('initcfg');
clear matlabbatch;

%file paths
base_filepath  = '/data/gpfs/projects/punim2118/bCFS/Data/';
spm_path = fullfile(base_filepath, 'preprocessed_files');
% output directory for DCM results
output_dir = fullfile(spm_path, 'dcm'); 
cd(output_dir);

p_names_all = {'S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S23', 'S24', 'S25', 'S26_b1','S26_b2', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33'};
% p_names is a clean variable that has S03, S23, and S26 b1 and b2 excluded
p_names = {'S01', 'S02', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S24', 'S25', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33'};
  

GCM = {};
for PP = 1:length(p_names)    
    dcm_filename = ['DCM_' p_names{PP} '.mat']; 
    dcm_path = fullfile(output_dir, dcm_filename);
    GCM = [GCM, dcm_path];
end
GCM = GCM';

gcm_filename = ['GCM.mat'];
fprintf('Saving fitted GCM: %s\n', gcm_filename);
save(fullfile(output_dir, gcm_filename), 'GCM');
