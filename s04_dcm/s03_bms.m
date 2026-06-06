%% bms family comparisons

% Kav Bandara  2025

clear; clc;

addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12');
spm('defaults', 'eeg');

% Base paths
base_dir = '/data/gpfs/projects/punim2118/bCFS/Data/preprocessed_files/dcm';
output_dir = fullfile(base_dir, 'bms_results');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

p_names_all = {'S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S23', 'S24', 'S25', 'S26_b1','S26_b2', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33'};

% p_names is a clean variable that has S03, S23, and S26 b1 and b2 excluded
p_names = {'S01', 'S02', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S24', 'S25', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33'};

n_subs = length(p_names);

models = {
    'M1',    'standard', fullfile(base_dir, 'M1', 'standard_model');
    'M1',    'ramping',  fullfile(base_dir, 'M1', 'ramping_model');
    'no_M1', 'standard', fullfile(base_dir, 'no_M1', 'standard_model');
    'no_M1', 'ramping',  fullfile(base_dir, 'no_M1', 'ramping_model');
};

n_models = size(models, 1);
model_names = {'M1_Standard', 'M1_Ramping', 'No_M1_Standard', 'No_M1_Ramping'};

subj = struct();

% setup model space 
for s = 1:n_subs
    for m = 1:n_models
        
        dcm_dir = models{m, 3};
        arch = models{m, 1};
        inp = models{m, 2};

        dcm_file = fullfile(dcm_dir, sprintf('DCM_%s_%s_%s.mat', arch, inp, p_names{s}));

        load(dcm_file, 'DCM');
        
        % subj(s).sess(1).model(m)
        subj(s).sess(1).model(m).fname = dcm_file;
        subj(s).sess(1).model(m).F = DCM.F;
        subj(s).sess(1).model(m).Ep = DCM.Ep;
        subj(s).sess(1).model(m).Cp = DCM.Cp;
        subj(s).sess(1).model(m).nonLin = 0;  
        
    end

end

% save 
model_space_file = fullfile(output_dir, 'model_space.mat');
save(model_space_file, 'subj');

%% create model families 

% standard vs ramping

family.names = {'Standard', 'Ramping'};
family.partition = [1, 2, 1, 2];

family_input_file = fullfile(output_dir, 'family_input_type.mat');
save(family_input_file, 'family');


% M1 vs no M1
clear family
family.names = {'M1', 'no_M1'};
family.partition = [1, 1, 2, 2];

family_arch_file = fullfile(output_dir, 'family_architecture.mat');
save(family_arch_file, 'family');

%% Run BMS on families 

% first BMS -- overall 

clear matlabbatch

matlabbatch = {};
matlabbatch{1}.spm.dcm.bms.inference.dir = {output_dir};
matlabbatch{1}.spm.dcm.bms.inference.sess_dcm = {};  
matlabbatch{1}.spm.dcm.bms.inference.model_sp = {model_space_file};
matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
matlabbatch{1}.spm.dcm.bms.inference.method = 'RFX';
matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};  
matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;  
matlabbatch{1}.spm.dcm.bms.inference.verify_id = 0;  

spm_jobman('run', matlabbatch);

% rename bms.mat to avoid overwriting
movefile(fullfile(output_dir, 'BMS.mat'), fullfile(output_dir, 'BMS_overall.mat'));

%% BMS 2 - compare input - standard vs ramping

clear matlabbatch

matlabbatch = {};
matlabbatch{1}.spm.dcm.bms.inference.dir = {output_dir};
matlabbatch{1}.spm.dcm.bms.inference.sess_dcm = {};
matlabbatch{1}.spm.dcm.bms.inference.model_sp = {model_space_file};
matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
matlabbatch{1}.spm.dcm.bms.inference.method = 'RFX';
matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {family_input_file};
matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
matlabbatch{1}.spm.dcm.bms.inference.verify_id = 0;

spm_jobman('run', matlabbatch);

% rename bms.mat to avoid overwriting
movefile(fullfile(output_dir, 'BMS.mat'), fullfile(output_dir, 'BMS_input_type.mat'));

%% BMS 3 - compare architecture - with M1 or without M1 

clear matlabbatch

matlabbatch = {};
matlabbatch{1}.spm.dcm.bms.inference.dir = {output_dir};
matlabbatch{1}.spm.dcm.bms.inference.sess_dcm = {};
matlabbatch{1}.spm.dcm.bms.inference.model_sp = {model_space_file};
matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
matlabbatch{1}.spm.dcm.bms.inference.method = 'RFX';
matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {family_arch_file};
matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
matlabbatch{1}.spm.dcm.bms.inference.verify_id = 0;

spm_jobman('run', matlabbatch);

% rename bms.mat to avoid overwriting
movefile(fullfile(output_dir, 'BMS.mat'), fullfile(output_dir, 'BMS_architecture.mat'));

