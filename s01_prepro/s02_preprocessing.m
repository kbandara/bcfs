%% Pre-Processing of EEG bCFS data 
% kav unimelb 2025

% this script takes data partially preprocessed from mne which has been 
% downsampled to 200 Hz, low pass filtered at 1hz and artefactual ica 
% components removed. then rereferenced to the average. 

% this script consists of 6 steps:
% 1) performs conversion into spm readable format
% 2) filtering and rereferencing
% 3) artefact detection - an intermediate step to load the correct triggers
%    from the raw files (details below)
% 4) epoching (this step includes a bunch of intermediate
%    steps too to clean the data including loading accuracy data from a
%    behavioural file)
% 5) marking data as bad and removing bad trials
% 6) robust averaging, a second low-pass filter, and then finally baseline 
%    correction. 

% this script allows for parallel processing (using parfor) as well -- just
% turn on at section 0  

clear all; close all

addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12')

spm('defaults', 'eeg');

%% %%%%%%%%           SECTION 0 - LOAD FILES & DEFINE VARIABLES          %%%%%%%%%%%%

% Settings and filenames for pre-processing
epochTimeWindow = [-100 3000]; %epoch around this time window (ms)

base_filepath  = '/data/gpfs/projects/punim2118/bCFS/Data/';
cd(base_filepath)

p_names_all = {'S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S23', 'S24', 'S25', 'S26_b1','S26_b2', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33'};

% p_names is a clean variable that has S03, S23, and S26 b1 and b2 excluded
p_names = {'S01', 'S02', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S24', 'S25', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33'};
%p_names = {'S01'};

new_folder = 'preprocessed_files';
if ~exist(fullfile(base_filepath, new_folder), 'dir')
    mkdir(base_filepath, new_folder)
end

%%              =====    PARALLEL PROCESSING SETTINGS    =====

% Set to true to use parallel processing, false to run serially
parallel_processing = true; 

% start parpool 
if parallel_processing
    % Check if a parallel pool is already running
    if isempty(gcp('nocreate'))
        disp('Attempting to start parallel pool...');
        
        % Check for Slurm-allocated CPUs
        n_cores_str = getenv('SLURM_CPUS_PER_TASK');
        
        if ~isempty(n_cores_str)
            % If running on Slurm, use the allocated cores
            n_cores = str2num(n_cores_str);
            fprintf('Slurm has allocated %d cores. Starting parpool with %d workers.\n', n_cores, n_cores);
            parpool(n_cores);
        else
            % If not on Slurm (e.g., running locally), start a default pool
            disp('Not running on Slurm. Starting a default parallel pool.');
            parpool;
        end
    else
        disp('Parallel pool already running.');
    end
end

%% %%%%%%%%%%             SECTION I - PREPROCESSING            %%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   STEP 1: Convert the datafile to SPM readable format   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize cell array for the jobs
convert_jobs = cell(1, length(p_names));

for c = 1:length(p_names)

    filename = fullfile(base_filepath, 'EEG', [p_names{c} '_bCFS_ica_filt.fif']); % 
      
    % Convert from EDF to dat/mat file
    matlabbatch = {};
    matlabbatch{1}.spm.meeg.convert.dataset = {filename};
    matlabbatch{1}.spm.meeg.convert.mode.continuous.readall = 1;
    matlabbatch{1}.spm.meeg.convert.channels{1}.all = 'all';
    matlabbatch{1}.spm.meeg.convert.outfile = fullfile(base_filepath, new_folder, ['spmeeg_' p_names{c} '_eeg.mat']);
    matlabbatch{1}.spm.meeg.convert.eventpadding = 0;
    matlabbatch{1}.spm.meeg.convert.blocksize = 3276800;
    matlabbatch{1}.spm.meeg.convert.checkboundary = 1;
    matlabbatch{1}.spm.meeg.convert.saveorigheader = 0;
    matlabbatch{1}.spm.meeg.convert.inputformat = 'autodetect';
    
    % save this participant's job to the main cell array
    convert_jobs{c} = matlabbatch;
end

%now run 
spm_jobman('initcfg');

if parallel_processing
   
    % Pre-allocate a cell array to store the outcome of each job
    job_outcomes = cell(1, length(p_names));

    parfor c = 1:length(p_names)
        try
            % Each worker runs one job from the pre-defined cell array
            spm_jobman('run', convert_jobs{c});
            job_outcomes{c} = 'Success';
        catch e
            % If a job fails, store 'failed' in the outcome cell
            job_outcomes{c} = sprintf('Failed: %s', e.message);
        end
    end
    
    % Display summary of outcomes after the parfor loop
    disp('--- Parallel processing complete. Job outcomes: ---');
    for c = 1:length(p_names)
        fprintf('Participant %s: %s\n', p_names{c}, job_outcomes{c});
    end

else
    spm_jobman('run', convert_jobs);
end

clear job_outcomes matlabbatch c

% now change path to where the new data file is 

for path = 1:length(p_names)
    
    fname = fullfile(base_filepath, new_folder, ['spmeeg_' p_names{path} '_eeg.mat']);
    
    load(fname, 'D');
    
    % update D.data.path
    [fpath, fbase, ~] = fileparts(fname);
    correct_dat_path = fullfile(fpath, [fbase '.dat']);
    
    % Fix the path
    D.data.fname = correct_dat_path;
    
    % Save the corrected struct back to file
    save(fname, 'D');
    
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%              STEP 2: Filter and Rereference             %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%high pass filter 

% Initialize cell array for the jobs
highpass_jobs = cell(1, length(p_names));

for fh = 1:length(p_names)

    filename = fullfile(base_filepath, new_folder, ['spmeeg_' p_names{fh} '_eeg.mat']);
    matlabbatch = {};
   
    matlabbatch{1}.spm.meeg.preproc.filter.D =  {filename};
    matlabbatch{1}.spm.meeg.preproc.filter.type = 'butterworth';
    matlabbatch{1}.spm.meeg.preproc.filter.band = 'high';
    matlabbatch{1}.spm.meeg.preproc.filter.freq = 0.1;
    matlabbatch{1}.spm.meeg.preproc.filter.dir = 'twopass';
    matlabbatch{1}.spm.meeg.preproc.filter.order = 5;
    matlabbatch{1}.spm.meeg.preproc.filter.prefix = 'fh_';
    
    % save this participant's job to the main cell array
    highpass_jobs{fh} = matlabbatch; 

end

spm_jobman('initcfg');

if parallel_processing
   
    % Pre-allocate a cell array to store the outcome of each job
    job_outcomes = cell(1, length(p_names));

    parfor i = 1:length(p_names)
        try
            % Each worker runs one job from the pre-defined cell array
            spm_jobman('run', highpass_jobs{i});
            job_outcomes{i} = 'Success';
        catch e
            % If a job fails, store 'failed' in the outcome cell
            job_outcomes{i} = sprintf('Failed: %s', e.message);
        end
    end
    
    % Display summary of outcomes after the parfor loop
    disp('--- Parallel processing complete. Job outcomes: ---');
    for i = 1:length(p_names)
        fprintf('Participant %s: %s\n', p_names{i}, job_outcomes{i});
    end

else
    spm_jobman('run', highpass_jobs);
end

clear job_outcomes matlabbatch i


% Remove line noise using a stopband filter at 50Hz (data collected in Brisbane, Australia)

% Initialize cell array for the jobs
linenoise_jobs = cell(1, length(p_names));

for ln = 1:length(p_names)

    filename = fullfile(base_filepath, new_folder, ['fh_spmeeg_' p_names{ln} '_eeg.mat']);
    matlabbatch = {};
    matlabbatch{1}.spm.meeg.preproc.filter.D = {filename};
    matlabbatch{1}.spm.meeg.preproc.filter.type = 'butterworth';
    matlabbatch{1}.spm.meeg.preproc.filter.band = 'stop';
    matlabbatch{1}.spm.meeg.preproc.filter.freq = [49 51]; % 50 Hz for Aus
    matlabbatch{1}.spm.meeg.preproc.filter.dir = 'twopass'; %default
    matlabbatch{1}.spm.meeg.preproc.filter.order = 5;
    matlabbatch{1}.spm.meeg.preproc.filter.prefix = 'ln_';
    
    % save this participant's job to the main cell array
    linenoise_jobs{ln} = matlabbatch; 
end

spm_jobman('initcfg');

if parallel_processing
   
    % Pre-allocate a cell array to store the outcome of each job
    job_outcomes = cell(1, length(p_names));

    parfor i = 1:length(p_names)
        try
            % Each worker runs one job from the pre-defined cell array
            spm_jobman('run', linenoise_jobs{i});
            job_outcomes{i} = 'Success';
        catch e
            % If a job fails, store 'failed' in the outcome cell
            job_outcomes{i} = sprintf('Failed: %s', e.message);
        end
    end
    
    % Display summary of outcomes after the parfor loop
    disp('--- Parallel processing complete. Job outcomes: ---');
    for i = 1:length(p_names)
        fprintf('Participant %s: %s\n', p_names{i}, job_outcomes{i});
    end

else
    spm_jobman('run', linenoise_jobs);
end

clear job_outcomes matlabbatch i

%% Remove 10Hz mask flicker using moving average filter

% moving average of 200ms to remove mask flicker. this is what was done in
% the original paper. note that these files were not used for decoding.

window_ms = 200; 

for ma_idx = 1:length(p_names)
    
    input_file = fullfile(base_filepath, new_folder, ['ln_fh_spmeeg_' p_names{ma_idx} '_eeg.mat']);
    
    D = spm_eeg_load(input_file);
    
    fs = D.fsample;
    window_samples = round(window_ms / 1000 * fs);
    
    if mod(window_samples, 2) == 0
        window_samples = window_samples + 1;
    end
    
    fprintf('%s: fs=%dHz, window=%d samples (%.1fms)\n', ...
            p_names{ma_idx}, fs, window_samples, window_samples/fs*1000);
    
    % data dimensions
    n_channels = D.nchannels;
    n_samples = D.nsamples;
    
    % raw data
    data = D(:,:,:); 
    
    if ndims(data) == 2
        smoothed_data = zeros(size(data));
        for ch = 1:n_channels
            smoothed_data(ch,:) = movmean(data(ch,:), window_samples);
        end
    end
    
    output_file = fullfile(base_filepath, new_folder, ['sm_ln_fh_spmeeg_' p_names{ma_idx} '_eeg.mat']);
    
    Dnew = clone(D, output_file, [n_channels, n_samples, size(data,3)]);
    Dnew(:,:,:) = smoothed_data;
    save(Dnew);
    
    clear D Dnew data smoothed_data
end

 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                 STEP 3: Artefact Detection               %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mark artefacts using a threshold of 100 uV

% Initialize cell array for the jobs
markartefacts_jobs = cell(1, length(p_names));

for ma = 1:length(p_names)

    filename = fullfile(base_filepath, new_folder, ['sm_ln_fh_spmeeg_' p_names{ma} '_eeg.mat']);
    matlabbatch = {};
    matlabbatch{1}.spm.meeg.preproc.artefact.D = {filename};
    matlabbatch{1}.spm.meeg.preproc.artefact.mode = 'mark';
    matlabbatch{1}.spm.meeg.preproc.artefact.badchanthresh = 0.2;
    matlabbatch{1}.spm.meeg.preproc.artefact.append = true;
    matlabbatch{1}.spm.meeg.preproc.artefact.methods(1).channels{1}.type = 'all';
    matlabbatch{1}.spm.meeg.preproc.artefact.methods(1).fun.threshchan.threshold = 100; %100 uV
    matlabbatch{1}.spm.meeg.preproc.artefact.methods(1).fun.threshchan.excwin = 1000;
    
    matlabbatch{1}.spm.meeg.preproc.artefact.prefix = 'a_';
    
    % save this participant's job to the main cell array
    markartefacts_jobs{ma} = matlabbatch; 
end

spm_jobman('initcfg');

if parallel_processing
   
    % Pre-allocate a cell array to store the outcome of each job
    job_outcomes = cell(1, length(p_names));

    parfor i = 1:length(p_names)
        try
            % Each worker runs one job from the pre-defined cell array
            spm_jobman('run', markartefacts_jobs{i});
            job_outcomes{i} = 'Success';
        catch e
            % If a job fails, store 'failed' in the outcome cell
            job_outcomes{i} = sprintf('Failed: %s', e.message);
        end
    end
    
    % Display summary of outcomes after the parfor loop
    disp('--- Parallel processing complete. Job outcomes: ---');
    for i = 1:length(p_names)
        fprintf('Participant %s: %s\n', p_names{i}, job_outcomes{i});
    end

else
    spm_jobman('run', markartefacts_jobs);
end

clear job_outcomes matlabbatch i

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%        INTERMEDIATE STEP: Replace D.trials struct        %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%this is a work around because the fif files from mne appear to have lost
%the events. but the raw edf files contain the right events when converted
%into an spm readable format. here the edf file is converted into spm, then
%the D.trials struct in the file we have been working with is replaced with
%D.trials struct from the raw file. then the intermediate file is deleted
%to conserve memory. 


% Initialize cell array for the jobs
raw_jobs = cell(1, length(p_names));

for raw = 1:length(p_names)

    filename = fullfile(base_filepath, 'EEG', [p_names{raw} '_bCFS.bdf']); 
      
    % Convert from EDF to dat/mat file
    matlabbatch = {};
    matlabbatch{1}.spm.meeg.convert.dataset = {filename};
    matlabbatch{1}.spm.meeg.convert.mode.continuous.readall = 1;
    matlabbatch{1}.spm.meeg.convert.channels{1}.all = 'all';
    matlabbatch{1}.spm.meeg.convert.outfile = fullfile(base_filepath, new_folder, ['spmeeg_' p_names{raw} '_raw.mat']);
    matlabbatch{1}.spm.meeg.convert.eventpadding = 0;
    matlabbatch{1}.spm.meeg.convert.blocksize = 3276800;
    matlabbatch{1}.spm.meeg.convert.checkboundary = 1;
    matlabbatch{1}.spm.meeg.convert.saveorigheader = 0;
    matlabbatch{1}.spm.meeg.convert.inputformat = 'autodetect';
    
    % save this participant's job to the main cell array
    raw_jobs{raw} = matlabbatch;
end

%now run 
spm_jobman('initcfg');

if parallel_processing
   
    % Pre-allocate a cell array to store the outcome of each job
    job_outcomes = cell(1, length(p_names));

    parfor c = 1:length(p_names)
        try
            % Each worker runs one job from the pre-defined cell array
            spm_jobman('run', raw_jobs{c});
            job_outcomes{c} = 'Success';
        catch e
            % If a job fails, store 'failed' in the outcome cell
            job_outcomes{c} = sprintf('Failed: %s', e.message);
        end
    end
    
    % Display summary of outcomes after the parfor loop
    disp('--- Parallel processing complete. Job outcomes: ---');
    for c = 1:length(p_names)
        fprintf('Participant %s: %s\n', p_names{c}, job_outcomes{c});
    end

else
    spm_jobman('run', raw_jobs);
end

%now that the files have been converted, take the D.trials struct and
%replace it with the faulty D.trials struct 

for fix = 1:length(p_names)

    raw_file = fullfile(base_filepath, new_folder, ['spmeeg_' p_names{fix} '_raw.mat']); 
    working_file = fullfile(base_filepath, new_folder, ['a_sm_ln_fh_spmeeg_' p_names{fix} '_eeg.mat']);

    load(raw_file); 
    raw = D.trials; 

    load(working_file); 

    D.trials = raw;

    save(working_file, "D")

    clear raw D
end

%now delete intermediate raw converted file to save space
   
disp('--- Deleting intermediate converted raw files... ---');

for p = 1:length(p_names)
    % Define the filenames for the intermediate .mat and .dat files
    mat_file = ['spmeeg_' p_names{p} '_raw.mat'];
    dat_file = ['spmeeg_' p_names{p} '_raw.dat'];
    
    % Construct the full file paths
    mat_filepath = fullfile(base_filepath, new_folder, mat_file);
    dat_filepath = fullfile(base_filepath, new_folder, dat_file);
    
    % Check if the .mat file exists and delete it
    if exist(mat_filepath, 'file')
        delete(mat_filepath);
        fprintf('Deleted: %s\n', mat_filepath);
    else
        fprintf('Warning: Could not find file to delete: %s\n', mat_filepath);
    end
    
    % Check if the .dat file exists and delete it
    if exist(dat_filepath, 'file')
        delete(dat_filepath);
        fprintf('Deleted: %s\n', dat_filepath);
    else
        fprintf('Warning: Could not find file to delete: %s\n', dat_filepath);
    end
end

clear job_outcomes matlabbatch c p


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                   STEP 4: Epoch trials                   %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
FYI: Trial Condition Indices​

1 = Expected Neutral  ​

2 = Unexpected Neutral​

3 = Expected Fearful  ​

4 = Unexpected Fearful
%}

% Initialize cell array for the jobs
epoch_jobs = cell(1, length(p_names));

for et = 1:length(p_names)

    filename = fullfile(base_filepath, new_folder, ['a_sm_ln_fh_spmeeg_' p_names{et} '_eeg.mat']);
    matlabbatch = {};

    matlabbatch{1}.spm.meeg.preproc.epoch.D = {filename};
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.timewin = epochTimeWindow;
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).conditionlabel = 'expected_neutral';
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventtype = 'STATUS';
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).eventvalue = 1; %expected neutral trigger values
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(1).trlshift = 0;
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).conditionlabel = 'unexpected_neutral';
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).eventtype = 'STATUS';
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).eventvalue = 2; %unexpected neutral trigger values
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(2).trlshift = 0;
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).conditionlabel = 'expected_fearful';
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).eventtype = 'STATUS';
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).eventvalue = 3; %expected fearful trigger values 
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(3).trlshift = 0;
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(4).conditionlabel = 'unexpected_fearful';
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(4).eventtype = 'STATUS';
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(4).eventvalue = 4; %unexpected fearful trigger values  
    matlabbatch{1}.spm.meeg.preproc.epoch.trialchoice.define.trialdef(4).trlshift = 0;
    matlabbatch{1}.spm.meeg.preproc.epoch.bc = 0;
    matlabbatch{1}.spm.meeg.preproc.epoch.eventpadding = 0;
    matlabbatch{1}.spm.meeg.preproc.epoch.prefix = 'e_';
 
    
    % save this participant's job to the main cell array
    epoch_jobs{et} = matlabbatch; 
end

spm_jobman('initcfg');

if parallel_processing
   
    % Pre-allocate a cell array to store the outcome of each job
    job_outcomes = cell(1, length(p_names));

    parfor i = 1:length(p_names)
        try
            % Each worker runs one job from the pre-defined cell array
            spm_jobman('run', epoch_jobs{i});
            job_outcomes{i} = 'Success';
        catch e
            % If a job fails, store 'failed' in the outcome cell
            job_outcomes{i} = sprintf('Failed: %s', e.message);
        end
    end
    
    % Display summary of outcomes after the parfor loop
    disp('--- Parallel processing complete. Job outcomes: ---');
    for i = 1:length(p_names)
        fprintf('Participant %s: %s\n', p_names{i}, job_outcomes{i});
    end

else
    spm_jobman('run', epoch_jobs);
end

clear job_outcomes matlabbatch i

%% Now take other behavioural measures from sXX_eegbehav.mat and map onto D.trials

% but first some cleaning to match the file structures 

% clean s01 because it seems some triggers were lost in the eeg file --

%{
%               WARNING!!! ONLY RUN THIS ONCE!!

load(fullfile(base_filepath, 'Behavioural' , 'S01_eegbehav.mat'));

%remove the first 70 trials 

eegbehav.noresponse = eegbehav.noresponse(78:end);
eegbehav.idx        = eegbehav.idx(78:end);
eegbehav.rt         = eegbehav.rt(78:end);
eegbehav.trials     = eegbehav.trials(78:end);
eegbehav.acc        = eegbehav.acc(78:end);

save(fullfile(base_filepath, 'Behavioural' , 'S01_eegbehav.mat'), 'eegbehav');
%}

%% now that behav files match continue with mapping accuracy and rt to eeg data

for acc = 1:length(p_names)

    filename = fullfile(base_filepath, new_folder, ['e_a_sm_ln_fh_spmeeg_' p_names{acc} '_eeg.mat']);
    load(filename)

    behav_file = fullfile(base_filepath, 'Behavioural' , [p_names{acc} '_eegbehav.mat']);
    load(behav_file)
    
    % write accuracy and RT data into D.trials as new field
    acc_cells = num2cell(eegbehav.acc);
    [D.trials(1:length(D.trials)).acc] = acc_cells{:};

    rt_cells = num2cell(eegbehav.rt);
    [D.trials(1:length(D.trials)).rt] = rt_cells{:};

    behav_trig = num2cell(eegbehav.trials);
    [D.trials(1:length(D.trials)).behav_trig] = behav_trig{:};
    
    save(filename, 'D');
    
    clear D
end

%% quality checks to make sure the mapping worked 

%check for proportion of matches by comparing D.trials.label with
%D.trials.behav_trig

% trigger to label mapping 
trig2lab = containers.Map( ...
    {1,                 2,                    3,                    4}, ...
    {'expected_neutral','unexpected_neutral','expected_fearful',   'unexpected_fearful'} );

%init
results = [];             
subj_tags = strings(0,1);    

for ii = 1:numel(p_names)                     
    subj_tag = p_names{ii};                    
    subj_num = sscanf(subj_tag,'S%f');        

    filename = fullfile(base_filepath, new_folder, ...
                        ['e_a_sm_ln_fh_spmeeg_' subj_tag '_eeg.mat']);
    if ~exist(filename,'file')
        results = [results; subj_num, 0, 0, 0, NaN];
        subj_tags(end+1,1) = string(subj_tag);
        continue
    end

    tmp = load(filename,'D');
    D   = tmp.D;

    nT = length(D.trials);
    behav_trig = nan(1,nT);
    labels     = strings(1,nT);

    for t = 1:nT
        % behav trig (1–4 or NaN)
        if isfield(D.trials(t),'behav_trig') && ~isempty(D.trials(t).behav_trig)
            behav_trig(t) = double(D.trials(t).behav_trig);
        else
            behav_trig(t) = NaN;
        end

        % label string 
        lab = "";
        if isfield(D.trials(t),'label') && ~isempty(D.trials(t).label)
            L = D.trials(t).label;
            if iscell(L),     L = L{1}; end
            lab = string(strtrim(L));
        end
        labels(t) = lab;
    end

    % expected label from behav trig
    expected = strings(1,nT);
    for t = 1:nT
        if ~isnan(behav_trig(t)) && isKey(trig2lab, behav_trig(t))
            expected(t) = string(trig2lab(behav_trig(t)));
        else
            expected(t) = string(missing);
        end
    end

    % valid comparisons (have both a behav trig and a label)
    valid   = ~ismissing(expected) & ~ismissing(labels);
    nValid  = sum(valid);

    match   = false(1,nT);
    match(valid) = strcmpi(labels(valid), expected(valid));
    nMatch  = sum(match);

    pctMatch = 100 * nMatch / max(nValid,1);

    results   = [results; subj_num, nT, nValid, nMatch, pctMatch];
    subj_tags(end+1,1) = string(subj_tag);

    fprintf('S%02d: trials=%d | valid=%d | matches=%d (%.1f%%)\n', ...
            subj_num, nT, nValid, nMatch, pctMatch);

    bad_ix = find(valid & ~match);
    if ~isempty(bad_ix)
        k = bad_ix(1:min(5,numel(bad_ix)));
        fprintf('  examples mismatched (trial: label vs expected):\n');
        for jj = 1:numel(k)
            fprintf('   - t=%d: "%s" vs "%s"\n', k(jj), labels(k(jj)), expected(k(jj)));
        end
    end
end

fprintf('\n=== SUMMARY ===\n');
valid_rows = ~isnan(results(:,5));
fprintf('Participants with valid compare: %d/%d\n', sum(valid_rows), size(results,1));
fprintf('Mean %% match: %.2f%% (SD %.2f)\n', mean(results(valid_rows,5)), std(results(valid_rows,5)));

%% DATA QUALITY CHECKS

% HERE TRIALS WITH INACCURATE RESPONSES AND RESPONSES WITHIN 500MS OF 
% STIMULUS PRESENTATION ARE MARKED AS BAD -- AND PROCESSED LATER 

rt_fast_thresh = 0.5;   % seconds (500 ms)

% --- for computing group means ---
nSub       = numel(p_names);
nT_all     = nan(nSub,1);
fast_cnt   = zeros(nSub,1);
noresp_cnt = zeros(nSub,1);
inc_cnt    = zeros(nSub,1);
bad_cnt    = zeros(nSub,1);

for qc = 1:length(p_names)

    fname = fullfile(base_filepath, new_folder, ['e_a_sm_ln_fh_spmeeg_' p_names{qc} '_eeg.mat']);
    load(fname, 'D');

    nT     = numel(D.trials);
    rt     = nan(nT,1);
    acc    = nan(nT,1);
    exbad  = false(nT,1);

    for t = 1:nT
        if isfield(D.trials(t),'rt')  && ~isempty(D.trials(t).rt),   rt(t)  = double(D.trials(t).rt);  end
        if isfield(D.trials(t),'acc') && ~isempty(D.trials(t).acc),  acc(t) = double(D.trials(t).acc); end
        if isfield(D.trials(t),'bad') && ~isempty(D.trials(t).bad),  exbad(t) = logical(D.trials(t).bad); end
    end

    % no response defined by rt being NaN
    noresp    = isnan(rt);

    % fast / incorrect dont apply to noresp trials
    fast      = (rt < rt_fast_thresh) & ~noresp;
    incorrect = (acc == 0)            & ~noresp;

    % existing OR fast OR incorrect
    newbad = exbad | fast | incorrect;

    for t = 1:nT
        D.trials(t).qc_fast      = fast(t);
        D.trials(t).qc_noresp    = noresp(t);
        D.trials(t).qc_incorrect = incorrect(t);

        if noresp(t)
            D.trials(t).label = 'no_response';
        end

        D.trials(t).bad = newbad(t);
    end

    % counts (for group means)
    nT_all(qc)     = nT;
    fast_cnt(qc)   = sum(fast);
    noresp_cnt(qc) = sum(noresp);
    inc_cnt(qc)    = sum(incorrect);
    bad_cnt(qc)    = sum(newbad);

    % overlap check (fast & incorrect can overlap)
    both_fi = sum(fast & incorrect);

    fprintf('%s: bad=%d/%d  (fast=%d, incorrect=%d; no_resp relabeled=%d; fast&incorrect=%d)\n', ...
        p_names{qc}, sum(newbad), nT, sum(fast), sum(incorrect), sum(noresp), both_fi);

    save(fname, 'D');
    clear D
end

% --- group means across participants ---
fprintf('\n=== QC MEANS ACROSS PARTICIPANTS ===\n');
fprintf('Mean counts: fast=%.1f, incorrect=%.1f, no_resp=%.1f, total_bad=%.1f\n', ...
    mean(fast_cnt,'omitnan'), mean(inc_cnt,'omitnan'), mean(noresp_cnt,'omitnan'), mean(bad_cnt,'omitnan'));

fprintf('Mean %% of trials: fast=%.2f%%, incorrect=%.2f%%, no_resp=%.2f%%, total_bad=%.2f%%\n', ...
    100*mean(fast_cnt./nT_all,'omitnan'), 100*mean(inc_cnt./nT_all,'omitnan'), ...
    100*mean(noresp_cnt./nT_all,'omitnan'), 100*mean(bad_cnt./nT_all,'omitnan'));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                     REMOVE BAD TRIALS                    %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(fullfile(base_filepath, new_folder));

for bads = 1:length(p_names)

    fname = fullfile(base_filepath, new_folder, ['e_a_sm_ln_fh_spmeeg_' p_names{bads} '_eeg.mat']); 
    load(fname, 'D');
    
    for t = 1:length(D.trials)
        if D.trials(t).qc_noresp == 1 || D.trials(t).rt > 2.7
            D.trials(t).bad = 1;
        end
    end

    save(fname, 'D');

    %remove these bad trials
    S = [];
    S.D = D.fname; 
    D = spm_eeg_remove_bad_trials(S);       

    clear D
end

% init cell array for the jobs
removebadtrials_jobs = cell(1, length(p_names));

for rr = 1:length(p_names)

    filename = fullfile(base_filepath, new_folder, ['re_a_sm_ln_fh_spmeeg_' p_names{rr} '_eeg.mat']);
    matlabbatch = {};
    matlabbatch{1}.spm.meeg.preproc.artefact.D = {filename};  
    matlabbatch{1}.spm.meeg.preproc.artefact.mode = 'reject';  
    matlabbatch{1}.spm.meeg.preproc.artefact.badchanthresh = 0.2; 
    matlabbatch{1}.spm.meeg.preproc.artefact.append = true;
    matlabbatch{1}.spm.meeg.preproc.artefact.methods(1).channels{1}.type = 'EEG';
    matlabbatch{1}.spm.meeg.preproc.artefact.methods(1).fun.events.whatevents.artefacts = 1; % 1 = "all"
    matlabbatch{1}.spm.meeg.preproc.artefact.prefix = 'r_';
    
    % save this participant's job to the main cell array
    removebadtrials_jobs{rr} = matlabbatch; 
end

spm_jobman('initcfg');

if parallel_processing
   
    % Pre-allocate a cell array to store the outcome of each job
    job_outcomes = cell(1, length(p_names));

    parfor i = 1:length(p_names)
        try
            % Each worker runs one job from the pre-defined cell array
            spm_jobman('run', removebadtrials_jobs{i});
            job_outcomes{i} = 'Success';
        catch e
            % If a job fails, store 'failed' in the outcome cell
            job_outcomes{i} = sprintf('Failed: %s', e.message);
        end
    end
    
    % Display summary of outcomes after the parfor loop
    disp('--- Parallel processing complete. Job outcomes: ---');
    for i = 1:length(p_names)
        fprintf('Participant %s: %s\n', p_names{i}, job_outcomes{i});
    end

else
    spm_jobman('run', removebadtrials_jobs);
end

clear job_outcomes matlabbatch i


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    STEP 6: Add breakthrough times from decoding to D.trials    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% to run this section the decoding scripts must have run already 

decoding_folder = fullfile(base_filepath, new_folder, 'decoding');

% load decoding results 
load(fullfile(decoding_folder, 'subject_voting_breakthrough.mat'), ...
     'subject_voting_breakthrough', 'rt_all', 'p_names');

for breakthru = 1:length(p_names)
    
    % Get this subject's voting breakthrough
    voting_bt = subject_voting_breakthrough(breakthru);
    mean_rt = rt_all(breakthru);
    
    % Load EEG file
    fname = fullfile(base_filepath, new_folder, ['r_re_a_sm_ln_fh_spmeeg_' p_names{breakthru} '_eeg.mat']);

    load(fname, 'D');
    
    % Add voting breakthrough to each trial (same value for all trials)
    for tr = 1:length(D.trials)
        D.trials(tr).voting_breakthrough = voting_bt;
        D.trials(tr).mean_rt = mean_rt;
    end
    
    % Store summary at file level for easy access
    D.breakthrough_summary.voting = voting_bt;
    D.breakthrough_summary.mean_rt = mean_rt;
    D.breakthrough_summary.bt_as_pct_rt = (voting_bt / mean_rt) * 100;
    
    save(fname, 'D');
    
    fprintf('%s: Voting BT = %.3fs (%.1f%% of RT)\n', ...
            p_names{breakthru}, voting_bt, (voting_bt / mean_rt) * 100);
    
    clear D
end

fprintf('\nSubject-level voting breakthrough added to all D.trials structs.\n');

%% assign channel coords - 3d position is missing so assign spm defaults

for coords = 1:length(p_names)

    fname = fullfile(base_filepath, new_folder, ['r_re_a_sm_ln_fh_spmeeg_' p_names{coords} '_eeg.mat']);
    
    D = spm_eeg_load(fname); 

    D = sensors(D, 'EEG', []); 
    
    eeg_indices = D.indchantype('EEG');
    n_eeg = length(eeg_indices);
    
    % create a matrix of NaNs (2 rows x n columns)
    D = coor2D(D, eeg_indices, NaN(2, n_eeg));
    
    save(fname, 'D');

    % assign new coords from library
    S = [];
    S.D = D;
    S.task = 'defaulteegsens'; 
    S.source = 'locs';         
    
    D = spm_eeg_prep(S);
    
    save(fname, 'D');

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 7: Robust average, low-pass filter and baseline correct  trials %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize cell array for the jobs
robustavg_jobs = cell(1, length(p_names));

for ra = 1:length(p_names)

    filename = fullfile(base_filepath, new_folder, ['r_re_a_sm_ln_fh_spmeeg_' p_names{ra} '_eeg.mat']);
    matlabbatch = {};
    matlabbatch{1}.spm.meeg.averaging.average.D = {filename};
    matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.ks = 3;
    matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.bycondition = true;
    matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.savew = false;
    matlabbatch{1}.spm.meeg.averaging.average.userobust.robust.removebad = true; 
    matlabbatch{1}.spm.meeg.averaging.average.plv = false;
    matlabbatch{1}.spm.meeg.averaging.average.prefix = 'ra_';
    
    % save this participant's job to the main cell array
    robustavg_jobs{ra} = matlabbatch; 
end

spm_jobman('initcfg');

if parallel_processing
   
    % Pre-allocate a cell array to store the outcome of each job
    job_outcomes = cell(1, length(p_names));

    parfor i = 1:length(p_names)
        try
            % Each worker runs one job from the pre-defined cell array
            spm_jobman('run', robustavg_jobs{i});
            job_outcomes{i} = 'Success';
        catch e
            % If a job fails, store 'failed' in the outcome cell
            job_outcomes{i} = sprintf('Failed: %s', e.message);
        end
    end
    
    % Display summary of outcomes after the parfor loop
    disp('--- Parallel processing complete. Job outcomes: ---');
    for i = 1:length(p_names)
        fprintf('Participant %s: %s\n', p_names{i}, job_outcomes{i});
    end

else
    spm_jobman('run', robustavg_jobs);
end

clear job_outcomes matlabbatch i

%low pass filter a second time (robust averaging introduces high frequency noise)
% Initialize cell array for the jobs
lowpass_jobs = cell(1, length(p_names));

for fl = 1:length(p_names)

    filename = fullfile(base_filepath, new_folder, ['ra_r_re_a_sm_ln_fh_spmeeg_' p_names{fl} '_eeg.mat']);
    matlabbatch = {};
    matlabbatch{1}.spm.meeg.preproc.filter.D = {filename};
    matlabbatch{1}.spm.meeg.preproc.filter.type = 'butterworth';
    matlabbatch{1}.spm.meeg.preproc.filter.band = 'low';
    matlabbatch{1}.spm.meeg.preproc.filter.freq = 40; % 40 Hz
    matlabbatch{1}.spm.meeg.preproc.filter.dir = 'twopass'; %default
    matlabbatch{1}.spm.meeg.preproc.filter.order = 5;
    matlabbatch{1}.spm.meeg.preproc.filter.prefix = 'fl_';
    
    % save this participant's job to the main cell array
    lowpass_jobs{fl} = matlabbatch; 
end

spm_jobman('initcfg');

if parallel_processing
   
    % Pre-allocate a cell array to store the outcome of each job
    job_outcomes = cell(1, length(p_names));

    parfor i = 1:length(p_names)
        try
            % Each worker runs one job from the pre-defined cell array
            spm_jobman('run', lowpass_jobs{i});
            job_outcomes{i} = 'Success';
        catch e
            % If a job fails, store 'failed' in the outcome cell
            job_outcomes{i} = sprintf('Failed: %s', e.message);
        end
    end
    
    % Display summary of outcomes after the parfor loop
    disp('--- Parallel processing complete. Job outcomes: ---');
    for i = 1:length(p_names)
        fprintf('Participant %s: %s\n', p_names{i}, job_outcomes{i});
    end

else
    spm_jobman('run', lowpass_jobs);
end

clear job_outcomes matlabbatch i

%baseline correction

% Initialize cell array for the jobs
baseline_jobs = cell(1, length(p_names));

for bc = 1:length(p_names)

    filename = fullfile(base_filepath, new_folder, ['fl_ra_r_re_a_sm_ln_fh_spmeeg_' p_names{bc} '_eeg.mat']);
    matlabbatch = {};
    matlabbatch{1}.spm.meeg.preproc.bc.D =  {filename};
    matlabbatch{1}.spm.meeg.preproc.bc.timewin = [-100 0];
    matlabbatch{1}.spm.meeg.preproc.bc.prefix = 'bc_';
    
    % save this participant's job to the main cell array
    baseline_jobs{bc} = matlabbatch; 
end

spm_jobman('initcfg');

if parallel_processing
   
    % Pre-allocate a cell array to store the outcome of each job
    job_outcomes = cell(1, length(p_names));

    parfor i = 1:length(p_names)
        try
            % Each worker runs one job from the pre-defined cell array
            spm_jobman('run', baseline_jobs{i});
            job_outcomes{i} = 'Success';
        catch e
            % If a job fails, store 'failed' in the outcome cell
            job_outcomes{i} = sprintf('Failed: %s', e.message);
        end
    end
    
    % Display summary of outcomes after the parfor loop
    disp('--- Parallel processing complete. Job outcomes: ---');
    for i = 1:length(p_names)
        fprintf('Participant %s: %s\n', p_names{i}, job_outcomes{i});
    end

else
    spm_jobman('run', baseline_jobs);
end

clear job_outcomes matlabbatch i


%%  Calculate aberage trials per participant
%before robust averaging
total_clean_trials = zeros(1, length(p_names));

for i = 1:length(p_names)
    
    filename = fullfile(base_filepath, new_folder, ['r_re_a_sm_ln_fh_spmeeg_' p_names{i} '_eeg.mat']);
    
    D = spm_eeg_load(filename);
    
    total_clean_trials(i) = D.ntrials;
    
end

% Calculate mean and standard deviation
mean_trials = mean(total_clean_trials);
sd_trials   = std(total_clean_trials);
