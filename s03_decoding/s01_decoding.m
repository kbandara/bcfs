%% s01_visual_decoding.m

% script to classify mask window and face window during bCFS 

% Kav Bandara, University of Melbourne, 2025


clear all; close all

addpath('/data/gpfs/projects/punim2118/envs/spm12/spm12')

spm('defaults', 'eeg');

%% config 

% random seed 
rng(42, 'twister');

base_filepath  = '/data/gpfs/projects/punim2118/bCFS/Data/preprocessed_files/';

cd(base_filepath)

output_folder = fullfile('decoding');

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

p_names_all = {'S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S23', 'S24', 'S25', 'S26_b1','S26_b2', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33'};

% p_names is a clean variable that has S03, S23, and S26 b1 and b2 excluded
p_names = {'S01', 'S02', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S24', 'S25', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33'};

% decoding params
mask_win = [0.2, 0.5];    % Class 0: mask  
face_win = [2.7, 3.0];    % Class 1: face  
search_win = [0.5, 2.7];

%cross validation 
n_folds = 7;             
kernel_scale = sqrt(64);
consecutive_samples = 5;
smooth_window = 10;      

%              =====    PARALLEL PROCESSING SETTINGS    =====

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

%% first some extra cleaning to get the files ready 

% we cant use the robustly averaged files so we will take intermediate
% files from preprocessing and apply baseline correction to them and also
% remove some extra trials (if responses happened after 2.7s)

% baseline correction

% Initialize cell array for the jobs
baseline_jobs = cell(1, length(p_names));

for bc = 1:length(p_names)

    filename = fullfile(base_filepath, ['r_e_a_ln_fh_spmeeg_' p_names{bc} '_eeg.mat']);
    matlabbatch = {};
    matlabbatch{1}.spm.meeg.preproc.bc.D =  {filename};
    matlabbatch{1}.spm.meeg.preproc.bc.timewin = [-100 0];
    matlabbatch{1}.spm.meeg.preproc.bc.prefix = 'bc_DECODING_FILES_';
    
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

clear matlabbatch i

%% mark trials with responses after 2.7s

% use my own logic because spm_remove_bad_trials removes all the info i
% have worked so hard to make - e.g. qc fields and rt fields

for bads = 1:length(p_names)

    fname = fullfile(base_filepath, ['bc_DECODING_FILES_r_e_a_ln_fh_spmeeg_' p_names{bads} '_eeg.mat']); 
    load(fname, 'D');
    
    for t = 1:length(D.trials)
        if D.trials(t).rt >= 2.7 || D.trials(t).qc_noresp == 1 
            D.trials(t).bad = 1;
        end
    end
    save(fname, 'D');

end

%% decoding 
parfor p = 1:length(p_names)
    
    sub_id = p_names{p};

    fname = fullfile(base_filepath, ['bc_DECODING_FILES_r_e_a_ln_fh_spmeeg_' sub_id '_eeg.mat']);
    S = load(fname); D = S.D;
    
    % Find EEG channel indices
    eeg_chan_idx = [];
    for ch = 1:length(D.channels)
        if strcmpi(D.channels(ch).type, 'EEG')
            eeg_chan_idx = [eeg_chan_idx, ch];
        end
    end

    % only EEG channels
    data = D.data(eeg_chan_idx, :, :);
    n_chans = length(eeg_chan_idx);
        
    % data dimensions
    n_samples = D.Nsamples;
    fs = D.Fsample;
    t_onset = D.timeOnset;
    
    time_axis = t_onset + (0:n_samples-1) / fs;
    
    [~, n_time, n_trials] = size(data);

    % true = good trial, false = bad trial
    good_trials = false(1, n_trials);
    for tr = 1:n_trials
        if isfield(D.trials(tr), 'bad')
            good_trials(tr) = (D.trials(tr).bad == 0);
        else
            good_trials(tr) = true;
        end
    end
    
    n_good = sum(good_trials);
    n_bad = sum(~good_trials);

    % use only good trials for training
    good_trial_idx = find(good_trials);
    
    %             --- CROSS VALIDATION ---

    % assign each trial to a block
    trials_per_block = n_trials / 14;
    block_indices = ceil((1:n_trials) / trials_per_block);
    block_indices = min(block_indices, 14); % avoids uneven divisibility issues 
    
    % pair blocks into folds (blocks 1-2 -> fold 1, blocks 3-4 -> fold 2)
    fold_indices = ceil(block_indices / 2);
    fold_indices = min(fold_indices, n_folds);
    
    % Convert time windows (in seconds) to sample indices
    
    % Mask window
    t_mask_idx = (time_axis >= mask_win(1)) & (time_axis <= mask_win(2));
    n_mask_samples = sum(t_mask_idx);
    
    % Face window
    t_face_idx = (time_axis >= face_win(1)) & (time_axis <= face_win(2));
    n_face_samples = sum(t_face_idx);
    
    % Search window
    t_search_idx = (time_axis >= search_win(1)) & (time_axis <= search_win(2));
       
    % matrix that stores predictions
    % note: each column is a trial, each row is a time point, and each value will be SVM scores (distance from decision boundary)
    score_matrix = NaN(n_time, n_trials);
    
    for fold = 1:n_folds
        fprintf('  Fold %d/%d: ', fold, n_folds);
        
        test_mask = (fold_indices == fold);   
        train_mask = ~test_mask & good_trials; 
        
        train_idx = find(train_mask);
        test_idx = find(test_mask);
        
        n_train = length(train_idx);
        n_test = length(test_idx);
        
        % training data
        X_train_raw = data(:, :, train_idx);  % channels x time x train_trials
        
        % get mask window time points
        X_0 = X_train_raw(:, t_mask_idx, :);  % channels x mask_samples x train_trials
        
        % Reshape so each row is one time point from one trial
        X_0 = permute(X_0, [2 3 1]);  
        
        % flatten the first two dimensions
        X_0 = reshape(X_0, [], n_chans);
        
        % label mask 
        Y_0 = zeros(size(X_0, 1), 1);  % Column of zeros
        
        % now repeat for faces
        X_1 = X_train_raw(:, t_face_idx, :);
        X_1 = permute(X_1, [2 3 1]);         
        X_1 = reshape(X_1, [], n_chans);   
        
        % label faces
        Y_1 = ones(size(X_1, 1), 1);  % Column of ones
        
        % --- Combine training data ---
        X_train = [X_0; X_1]; 
        Y_train = [Y_0; Y_1];
        
        % now standardise (z-score)
        % z-score for train data, then apply transformation to test data
        
        %mean and sd 
        mu = mean(X_train, 1); 
        sigma = std(X_train, 0, 1); 
        sigma(sigma == 0) = 1;
        
        % now z score 
        X_train = (X_train - mu) ./ sigma;
        
        % --- NOW Train SVM ---
                
        SVMModel = fitcsvm(X_train, Y_train, ...
            'KernelFunction', 'rbf', ...
            'KernelScale', kernel_scale, ...
            'ClassNames', [0, 1], ...
            'Standardize', false, ...
            'Prior', 'uniform');

        % CHANGED TO LINEAR 
        %SVMModel = fitcsvm(X_train, Y_train, ...
        %    'KernelFunction', 'linear', ...  
        %    'ClassNames', [0, 1], ...
        %    'Standardize', false, ...         
        %    'Prior', 'uniform');

        % NOW TEST!

        % Get test data (ALL time pointsm and trials including bad ones)
        X_test_raw = data(:, :, test_idx);          
        [~, n_time_local, n_test_trials] = size(X_test_raw);
        
        % reshape
        X_test_perm = permute(X_test_raw, [2 3 1]);  
        X_test_flat = reshape(X_test_perm, [], n_chans); 
        
        % apply same standardisation
        X_test_flat = (X_test_flat - mu) ./ sigma;
        
        % Get predictions
        % IF LINEAR: face_scores_flat = X_test_flat * SVMModel.Beta + SVMModel.Bias;
        [~, face_scores_raw] = predict(SVMModel, X_test_flat);
        face_scores_flat = face_scores_raw(:, 2);  
        
        face_scores = reshape(face_scores_flat, n_time_local, n_test_trials);
        
        % Store in the master matrix
        score_matrix(:, test_idx) = face_scores;

    end

    % some post processing 

    %baseline correcting 
    baseline_idx = (time_axis >= mask_win(1)) & (time_axis <= mask_win(2));
    trial_baselines = nanmean(score_matrix(baseline_idx, :), 1); 
    corrected_scores = score_matrix - trial_baselines;
    
    % calculate similarity 
    mean_score = nanmean(corrected_scores(:, good_trials), 2);

    mask_reference = nanmean(mean_score(t_mask_idx));   
    face_reference = nanmean(mean_score(t_face_idx)); 
    
    % scale to 0-100%
    score_range = face_reference - mask_reference;
    if score_range > 0
        face_similarity = ((mean_score - mask_reference) / score_range) * 100;
    else
        face_similarity = 50 * ones(size(mean_score));
    end
    
    face_similarity = max(min(face_similarity, 120), -20);
    
    smoothed_scores = smoothdata(corrected_scores, 1, 'gaussian', smooth_window);
    
    % compute bt times on smoothed scores 
    breakthrough_times = NaN(n_trials, 1);
    search_time_indices = find(t_search_idx);
    
    for tr = 1:n_trials
        if ~good_trials(tr)
            continue;
        end
        
        curve = smoothed_scores(t_search_idx, tr);
        
        above_threshold = curve > 0; % when curve crosses threshold
              
        kernel = ones(consecutive_samples, 1);
        run_sums = conv(double(above_threshold), kernel, 'valid');
        stable_start = find(run_sums == consecutive_samples, 1, 'first'); % 5 samples above threshold
        
        if ~isempty(stable_start)

            orig_idx = search_time_indices(stable_start);
            breakthrough_times(tr) = time_axis(orig_idx);
        end
    end

    % SAVE 
    save_name = fullfile(output_folder, [sub_id '_visual_decoding_results.mat']);
    parsave(save_name, score_matrix, corrected_scores, smoothed_scores, face_similarity, ...
     breakthrough_times, time_axis, good_trials, ...
     mask_win, face_win, search_win, consecutive_samples, ...
     mask_reference, face_reference, score_range);fprintf('  Saved: %s\n', save_name);
    
end

% save function for parfor 
function parsave(save_name, score_matrix, corrected_scores, smoothed_scores, face_similarity, ...
     breakthrough_times, time_axis, good_trials, ...
     mask_win, face_win, search_win, consecutive_samples, ...
     mask_reference, face_reference, score_range)
    save(save_name, 'score_matrix', 'corrected_scores', 'smoothed_scores', 'face_similarity', ...
     'breakthrough_times', 'time_axis', 'good_trials', ...
     'mask_win', 'face_win', 'search_win', 'consecutive_samples', ...
     'mask_reference', 'face_reference', 'score_range', '-v7.3');
end