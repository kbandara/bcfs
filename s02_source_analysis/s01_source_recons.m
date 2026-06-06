%% Source inversion for EEG/bCFS data
% Kav Bandara, University of Melbourne, 2025

%{

this script runs source inversion using settings recommended by SPM. Using
SPM batch scripting, this code involves: 
   1) source space modelling
   2) data co-registration
   3) forward modelling
   4) inverse reconstruction
   5) create source reconstructed images

%}

clear all; close all

addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12')

spm('defaults', 'eeg');

%%%%%%%%%%           SECTION 0 - LOAD FILES & DEFINE VARIABLES          %%%%%%%%%%%%

%% Settings and filenames for pre-processing

% Settings and filenames for pre-processing
epochTimeWindow = [500 1250]; %time window of interest

base_filepath  = '/data/gpfs/projects/punim2118/bCFS/Data/';

preprocessed_path = fullfile(base_filepath, 'preprocessed_files');
cd(preprocessed_path)

image_dir = fullfile(base_filepath, 'source_contrasts');

prefix = 'bc_fl_ra_r_re_a_sm_ln_fh_spmeeg_';
          
% Create the centralized output directory if it doesn't exist
if ~exist(image_dir, 'dir')
    fprintf('Creating directory: %s\n', image_dir);
    mkdir(image_dir);
end

cd(image_dir)

new_folder = 'group_GLM';
GLM_dir = fullfile(image_dir, new_folder);  

p_names_all = {'S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S23', 'S24', 'S25', 'S26_b1','S26_b2', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33'};

% p_names is a clean variable that has S03, S23, and S26 b1 and b2 excluded
p_names = {'S01', 'S02', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S24', 'S25', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33'};
%p_names = {'S01'};

roi_labels = {'L_V1', 'R_V1', 'L_FFA', 'R_FFA', 'L_M1', 'R_M1','L_PFC', 'R_PFC'};

%% fix missing elec positions
% Loop through your participants
for coords = 1:length(p_names)

    fname = fullfile(preprocessed_path, [ prefix, p_names{coords} '_eeg.mat']);
    
    D = spm_eeg_load(fname); 

    % Wipe 3D Sensors
    D = sensors(D, 'EEG', []); 
    
    % Wipe 2D Layout 
    eeg_indices = D.indchantype('EEG');
    n_eeg = length(eeg_indices);
    D = coor2D(D, eeg_indices, NaN(2, n_eeg));
    
    % Assign New Standard Coordinates
    S = [];
    S.D = D;
    S.task = 'defaulteegsens'; % Assign default 3D positions
    S.source = 'locs';         % Use standard library locations
    D = spm_eeg_prep(S);

    % Manually add standard fiducials
    % These are standard MNI coordinates for fiducials
    fid = [];
    fid.fid.label = {'Nasion', 'LPA', 'RPA'}';
    fid.fid.pnt = [
        0,    84,  -28;    % Nasion
       -83,  -20,  -65;    % LPA (left preauricular)
        83,  -20,  -65     % RPA (right preauricular)
    ];
    fid.unit = 'mm';
    fid.pnt = [];  % No headshape points
    
    D = fiducials(D, fid);
    % Final Save
    save(fname, 'D');

end

%% RUN SOURCE INVERSION - GROUP INVERSION 

group_inv = cell(length(p_names), 1);
for i = 1:length(p_names)
    group_inv{i} = fullfile(preprocessed_path, [prefix, p_names{i} '_eeg.mat']);
end

matlabbatch = {};
%Source space modelling (using templates), Coregister, Forward Model
matlabbatch{1}.spm.meeg.source.headmodel.D = group_inv;
matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1; %create mesh using sMRI 
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;    
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'LPA';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'lpa'; % enter 1 × 3 vector of 3d coordinates for each fiducial
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'Nasion';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'nas'; % enter 1 × 3 vector of 3d coordinates for each fiducial
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'RPA';
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa'; % enter 1 × 3 vector of 3d coordinates for each fiducial
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 0; % 0 = 'no, no is advised for MEG with sMRI; 1 is advised for EEG
matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';

%Model inversion and specify time (and frequency) window, using Multiple Sparse Priors algorithm 
matlabbatch{2}.spm.meeg.source.invert.D = group_inv;
matlabbatch{2}.spm.meeg.source.invert.val = 1; %save inversion to new index
matlabbatch{2}.spm.meeg.source.invert.whatconditions.all = 1; %use all conditions
matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.invtype = 'GS'; %Minimum norm is IID %Multiple sparse priors (greedy search)
matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.woi = epochTimeWindow; %time window of interest 
matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.foi = [0 256];%;[0 256]; %frequency window of interest 
matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.hanning = 0; %1 = Hanning taper at start and end of trial (0 = no taper, 1 = yes)
matlabbatch{2}.spm.meeg.source.invert.modality = {'EEG'};

%Display results within time (and frequency) window and create images
matlabbatch{3}.spm.meeg.source.results.D = group_inv;
matlabbatch{3}.spm.meeg.source.results.val = 1; %display inversion results from index
matlabbatch{3}.spm.meeg.source.results.woi = epochTimeWindow; % time of interest
matlabbatch{3}.spm.meeg.source.results.foi = [0 256]; % frequency window specify
matlabbatch{3}.spm.meeg.source.results.ctype = 'evoked'; % evoked is the option for averaged trials
matlabbatch{3}.spm.meeg.source.results.space = 1; % 1 = MNI or Native
matlabbatch{3}.spm.meeg.source.results.format = 'image';
matlabbatch{3}.spm.meeg.source.results.smoothing = 12; %Smoothing mm^3

spm_jobman('run', matlabbatch);
clear matlabbatch;   

%% MOVE FILES TO ONSET CONTRASTS FOLDER 
cd(preprocessed_path)
for move = 1:length(p_names)

    for i = 1:4 %for each condition/image generated above 
    
    og_file = fullfile(preprocessed_path, [prefix, p_names{move} '_eeg_1_t500_1250_f0_256_' num2str(i) '.nii']);
    target_file = fullfile(image_dir, [prefix, p_names{move} '_eeg_1_t500_1250_f0_256_' num2str(i) '.nii']);
       
    %move the files bc_fl_ra_r_e_a_ln_fh_spmeeg_S01
    movefile(og_file, target_file, 'f'); % 'f' to force overwrite if target exists

    clear og_file target_file; 

    end
end

%% RENAME SOURCE RECONSTRUCTED IMAGES WITH CONDITION LABELS from D.Trials SCRIPT

%%%%%%%%%%           SECTION 0 - DEFINE VARIABLES          %%%%%%%%%%%%

% names to construct each participants mat file and nii image file paths
nii_suffix = '_eeg_1_t500_1250_f0_256_';

cd(image_dir)

for i = 1:numel(p_names)
    
    current_subject = p_names{i};
    fprintf('\nProcessing Subject: %s\n', current_subject);
    
    eeg_filename = [prefix, current_subject, '_eeg.mat'];
    eeg_filepath = fullfile(preprocessed_path, eeg_filename);
 
    load(eeg_filepath);
    
    % Get the list of condition labels
    condition_labels = {D.trials.label};
    fprintf('Found %d conditions for this subject.\n', numel(condition_labels));

    % --- Loop through conditions and rename files ---
    for j = 1:numel(condition_labels)
        
        label = condition_labels{j};
        
        % Construct the old and new filenames
        old_filename = sprintf('%s%s%s%d.nii', prefix, current_subject, nii_suffix, j);
        old_filepath = fullfile(image_dir, old_filename);

        new_filename = [current_subject, '_500-1250ms_', label, '.nii'];
        new_filepath = fullfile(image_dir, new_filename);

        fprintf(' -> Renaming: %s  TO  %s\n', old_filename, new_filename);
        copyfile(old_filepath, new_filepath);
        
    end
end

%% RUN GLM CONTRASTS OF SOURCE RECONSTRUCTED IMAGES 

unexpected = '_unexpected';
expected   = '_expected';
fearful    = '_fearful';
neutral    = '_neutral';

new_folder = 'group_GLM';
GLM_dir = fullfile(image_dir, new_folder);  

if exist(GLM_dir, 'dir'); rmdir(GLM_dir, 's'); end
mkdir(GLM_dir);

UF = {}; % Unexpected Fearful
UN = {}; % Unexpected Neutral
EF = {}; % Expected Fearful
EN = {}; % Expected Neutral

% get images
for i = 1:length(p_names)
    
    sub = p_names{i};
    
    pattern = fullfile(image_dir, [sub, '*.nii']);
    d = dir(pattern);
    all_files = fullfile({d.folder}, {d.name})'; % Get full paths
       
    % logical indexing to find specific conditions
    is_unexp = contains(all_files, unexpected);
    is_exp   = contains(all_files, expected);
    is_fear  = contains(all_files, fearful);
    is_neut  = contains(all_files, neutral);
    
    % Sort into lists
    UF = [UF; all_files(is_unexp & is_fear)];
    UN = [UN; all_files(is_unexp & is_neut)];
    EF = [EF; all_files(is_exp   & is_fear)];
    EN = [EN; all_files(is_exp   & is_neut)];
end

clear matlabbatch;

matlabbatch{1}.spm.stats.factorial_design.dir = {GLM_dir};

matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Expectation';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1; % 1=Dependent (Within-subject)
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 0; 
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'Emotion';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1; % 1=Dependent
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 0; 
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;

% Cell 1,1: Unexpected Fearful
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1 1]; 
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans = UF;

% Cell 1,2: Unexpected Neutral
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [1 2]; 
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = UN;

% Cell 2,1: Expected Fearful
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [2 1]; 
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = EF;

% Cell 2,2: Expected Neutral
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [2 2]; 
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans = EN;

matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spmmat_file = fullfile(GLM_dir, 'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.spmmat = {spmmat_file};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat = {spmmat_file};
matlabbatch{3}.spm.stats.con.delete = 1;

% main effect of mean (All > Baseline)
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'All > Baseline';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 1 1 1]; 
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

% main effect of emotion (Fearful > Neutral)
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Fear > Neutral';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 -1 1 -1]; 
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

% main effect of expectation (Unexpected > Expected)
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'Unexpected > Expected';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [1 1 -1 -1];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

% interaction
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'Interaction';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [1 -1 -1 1];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

spm_jobman('run', matlabbatch);

%% CREATE RESULTS! EXTRACT PEAK MNI COORDINATES FROM GLM CONTRASTS

clear matlabbatch

outfile = fullfile(GLM_dir, 'indiv_mni_coords.m');

% ROI labels in order 
roi_labels = {'L_V1', 'R_V1', 'L_FFA', 'R_FFA', 'L_M1', 'R_M1','L_PFC', 'R_PFC'};

spmmat_file = fullfile(GLM_dir, 'SPM.mat'); 
    
for i = 1:length(roi_labels)

    %path to mask for this roi
    mask_dir = '/data/gpfs/projects/punim2118/bCFS/roi_masks/jubrain_masks'; %need to create these masks from jubrain (anatomy toolbox) first 
    mask_file = fullfile(mask_dir, [roi_labels{i} '.nii']);

    matlabbatch = {};

    %create results 
    matlabbatch{1}.spm.stats.results.spmmat = {spmmat_file};
    matlabbatch{1}.spm.stats.results.conspec.titlestr = roi_labels{i};
    matlabbatch{1}.spm.stats.results.conspec.contrasts = 1; % specify the contrast index -- im pretty sure this is 2? see SPM.xCon 
    matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none'; %none is ok because multiple comparisons isn't a problem here - we just want to find where there is the most robust activation
    matlabbatch{1}.spm.stats.results.conspec.thresh = 1.0; %note 0.05 is 95%. setting to 1 means spm will get any voxel that is different/positive for a participant 
    matlabbatch{1}.spm.stats.results.conspec.extent = 0;
    matlabbatch{1}.spm.stats.results.conspec.conjunction = 1; 
    %matlabbatch{1}.spm.stats.results.conspec.mask.none = 1; %if no mask 
    matlabbatch{1}.spm.stats.results.conspec.mask.image.name = {mask_file}; %filepath to mask .nii
    matlabbatch{1}.spm.stats.results.conspec.mask.image.mtype = 0; %0 is inclusive mask 
    matlabbatch{1}.spm.stats.results.units = 1;
    matlabbatch{1}.spm.stats.results.export{1}.ps = true;
    matlabbatch{1}.spm.stats.results.export{2}.jpg = true;
    matlabbatch{1}.spm.stats.results.export{3}.csv = true;
    
    startTime = now; % get a timestamp -- used for renaming files later

    spm_jobman('run', matlabbatch);

    % Find all CSV files matching the SPM default format in the GLM directory
    spm_csv_files = dir(fullfile(GLM_dir, 'spm_*.csv'));
    [~, latest_idx] = max([spm_csv_files.datenum]); %get latest file        
    old_filename = fullfile(GLM_dir, spm_csv_files(latest_idx).name);
    new_filename = fullfile(GLM_dir, [roi_labels{i} '.csv']);
    movefile(old_filename, new_filename); %now rename        
    fprintf('Renamed %s to %s\n', spm_csv_files(latest_idx).name, [roi_labels{i} '.csv']);

    clear matlabbatch

end

%%                    --- Get Peak MNI Coords ---

% init
indiv_mni_coords = cell(length(p_names), 2);

% create a 3xN matrix to hold the 3D coordinates for this participant
roi_coords_matrix = nan(3, length(roi_labels));

% loop through  ROIs
for r = 1:length(roi_labels)
    
    this_roi = roi_labels{r};
    csv_file = fullfile(GLM_dir, [this_roi '.csv']);

    % read the data from the CSV, skipping the 2 header lines.
    data = readmatrix(csv_file, 'HeaderLines', 2);

    % The x,y,z coordinates are in columns 12, 13, and 14.
    % We only want the first row - i.e. the peak.
    peak_coord = data(1, 12:14)'; 
    roi_coords_matrix(:, r) = peak_coord;
end 

indiv_mni_coords{1, 1} = roi_coords_matrix;

% Save the final cell array to a .mat file

output_mat_file = fullfile(GLM_dir, 'indiv_peak_coords.mat');

save(output_mat_file, 'indiv_mni_coords');
