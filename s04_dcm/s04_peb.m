%% s04_peb

% script to run peb analysis on bCFS DCMs

% kav bandara unimelb 2025

clear; close all; clc;

%% setup

addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12');
spm('defaults', 'eeg');

spmeeg_path = '/data/gpfs/projects/punim2118/bCFS/Data/preprocessed_files';
decoding_path = fullfile(spmeeg_path, 'decoding');
dcm_path = fullfile(spmeeg_path, 'dcm', 'M1', 'ramping_model'); %change to winning DCM model -- M1 ramping DCMs
peb_path = fullfile(spmeeg_path, 'dcm', 'peb'); if ~exist(peb_path, 'dir'); mkdir(peb_path); end
cd(dcm_path);

p_names_all = {'S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S23', 'S24', 'S25', 'S26_b1','S26_b2', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33'};
% p_names is a clean variable that has S03, S23, and S26 b1 and b2 excluded
p_names = {'S01', 'S02', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S24', 'S25', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33'};

n_subs = length(p_names);

model_to_test = 'M1_ramping';

condition_names = {'exp_neu', 'unexp_neu', 'exp_fear', 'unexp_fear'};

%% create GCM

GCM = {};
for PP = 1:length(p_names)    
    dcm_filename = ['DCM_' model_to_test '_' p_names{PP} '.mat'];
    GCM = [GCM, fullfile(dcm_path, dcm_filename)];
end
GCM = GCM';

% save
gcm_filename = ['GCM_' model_to_test '.mat'];
save(fullfile(peb_path, gcm_filename), 'GCM');

%% run PEB on A and B matrices

gcm_file = fullfile(peb_path, ['GCM_' model_to_test '.mat']);
load(gcm_file, 'GCM');

% Run on A-matrix 
M = struct();
M.Q = 'all';
M.X = ones(n_subs, 1);
field = {'A'};

PEB = spm_dcm_peb(GCM, M, field);
BMA = spm_dcm_peb_bmc(PEB);

save(fullfile(peb_path, ['PEB_A_' model_to_test '.mat']), 'PEB');
save(fullfile(peb_path, ['BMA_A_' model_to_test '.mat']), 'BMA');

clear BMA PEB M field

% Run on B-matrix

M = struct();
M.Q = 'all';
M.X = ones(n_subs, 1);
field = {'B{2}'};

PEB = spm_dcm_peb(GCM, M, field);
BMA = spm_dcm_peb_bmc(PEB);

save(fullfile(peb_path, ['PEB_B_' model_to_test '.mat']), 'PEB');
save(fullfile(peb_path, ['BMA_B_' model_to_test '.mat']), 'BMA');

clear BMA PEB M field

%% Run nested peb on A matrix PFC  

load(fullfile(peb_path, ['PEB_A_' model_to_test '.mat']), 'PEB');
load(fullfile(peb_path, ['GCM_' model_to_test '.mat']), 'GCM');

%get template from a DCM -- template is the full model
load_DCM = GCM{1};
load(load_DCM);
DCM_full = DCM;
if isfield(DCM_full, 'M')
    DCM_full = rmfield(DCM_full, 'M');
end
   
DCM_PFC = DCM_full;

DCM_PFC.A{2}(3,7) = 0; % L_FG <- L_PFC
DCM_PFC.A{2}(4,8) = 0; % R_FG <- R_PFC
DCM_PFC.A{2}(5,7) = 0; % L_M1 <- L_PFC
DCM_PFC.A{2}(6,8) = 0; % R_M1 <- R_PFC
DCM_PFC.A{2}(2,8) = 0; % R_V1 <- R_PFC
DCM_PFC.A{2}(1,7) = 0; % L_V1 <- L_PFC

% Run BMC
[BMA_PFC, BMR_PFC] = spm_dcm_peb_bmc(PEB, {DCM_full, DCM_PFC});

% Extract results
F_full = BMR_PFC{1}.F;
F_reduced = BMR_PFC{2}.F;
dF = F_full - F_reduced;

PFC_df = dF;
PFC_prob = 1 / (1 + exp(-dF));

fprintf('No feedback PFC: dF = %.2f, Pp = %.3f\n', dF, PFC_prob);

clear GCM PEB

% Save results
save(fullfile(peb_path, 'PFC_A_reduced_PEB.mat'), 'PFC_prob', 'PFC_df');

%% RUN NESTED PEB ON EXPECTATION EFFECT 

load(fullfile(peb_path, ['PEB_B_' model_to_test '.mat']), 'PEB');
load(fullfile(peb_path, ['GCM_' model_to_test '.mat']), 'GCM');

%get template from a DCM -- template is the full model
load_DCM = GCM{1};
load(load_DCM);
DCM_full = DCM;
if isfield(DCM_full, 'M')
    DCM_full = rmfield(DCM_full, 'M');
end
   
% No-PFC intrinsic model (remove intrinsic PFC connections)
DCM_PFC = DCM_full;

DCM_PFC.B{2}(3,7) = 0; % L_FG <- L_PFC
DCM_PFC.B{2}(4,8) = 0; % R_FG <- R_PFC
DCM_PFC.B{2}(5,7) = 0; % L_M1 <- L_PFC
DCM_PFC.B{2}(6,8) = 0; % R_M1 <- R_PFC
DCM_PFC.B{2}(2,8) = 0; % R_V1 <- R_PFC
DCM_PFC.B{2}(1,7) = 0; % L_V1 <- L_PFC

% Run BMC
[BMA_PFC, BMR_PFC] = spm_dcm_peb_bmc(PEB, {DCM_full, DCM_PFC});

% Extract results
F_full = BMR_PFC{1}.F;
F_reduced = BMR_PFC{2}.F;
dF = F_full - F_reduced;

PFC_df = dF;
PFC_prob = 1 / (1 + exp(-dF));

fprintf('No feedback PFC: dF = %.2f, Pp = %.3f\n', dF, PFC_prob);

clear GCM PEB

% Save results
save(fullfile(peb_path, 'PFC_B_expectation_PEB.mat'), 'PFC_prob', 'PFC_df');

%% REVIEW: PEB results in GUI (interactive step)

clear PEB GCM BMA 
%load GCM & PEB you want to view results for 

load(fullfile(peb_path, ['PEB_B_' model_to_test '.mat']), 'PEB');
load(fullfile(peb_path, ['GCM_' model_to_test '.mat']), 'GCM');

% Search over nested PEB models.
BMA = spm_dcm_peb_bmc(PEB);

spm_dcm_peb_review(BMA,GCM)
