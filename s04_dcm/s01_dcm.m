function s01_dcm(p_name, input_type)
    %{
    % this takes a time window from the SLURM scheduler to run in multiple
    % job arrays - i.e. parallel - with one participant per node

    % here the function takes arguments from the slurm job specifying the
    % p_name and input_type (standard/ramping)

    kav 2025
    unimelb
    
    %}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                                                         %
%                                                                         %
%                            TESTING BLOCK                                
%                                                                         %
%                                                                         %
% input_type = 'standard'; p_name = 'S01';  
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(input_type, 'standard')
        spm_root = '/data/gpfs/projects/punim2118/envs/spm12_standard/spm12';
    else
        spm_root = '/data/gpfs/projects/punim2118/envs/spm12_ramping/spm12';
    end

    addpath(spm_root);
    spm('defaults', 'eeg');
    spm_jobman('initcfg');

    base_filepath  = '/data/gpfs/projects/punim2118/bCFS/Data/';

    % dir containing spmeeg files 
    spmeeg_path = fullfile(base_filepath, 'preprocessed_files');
    cd(spmeeg_path);

    % Location of custom spm_erp_u files -- note this is not used in this
    % version because we replace the spm folder itself
    custom_functions_dir = '/data/gpfs/projects/punim2118/bCFS/Scripts/matlab/dcm/custom_input_functions/';

    % participant names - for debugging. when running as slurm job, participants are inputted from the slurm job via the bash script
    % p_names = {'S01', 'S02', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S24', 'S25', 'S27', 'S28', 'S29', 'S30', 'S31', 'S32', 'S33'};
      
    %load each participants peak mni coordinates 
    indiv_mni_coords = fullfile(base_filepath, 'source_contrasts', 'group_GLM', 'indiv_peak_coords.mat');
    load(indiv_mni_coords);
    indiv_mni_coords = indiv_mni_coords{1,1};

    %load decoding results
    decoding_results_file = fullfile(base_filepath, 'preprocessed_files', 'decoding', 'subject_voting_breakthrough.mat');
    decoding_data = load(decoding_results_file, 'subject_voting_breakthrough', 'p_names');
    p_idx = find(strcmp(decoding_data.p_names, p_name));
    bt_time_sec = decoding_data.subject_voting_breakthrough(p_idx);
    dcm_end_time = round(bt_time_sec * 1000); %convert to ms
    % ROI labels in order 
    roi_labels = {'L_V1', 'R_V1', 'L_FFA', 'R_FFA', 'L_M1', 'R_M1', 'L_PFC', 'R_PFC'};

    % output directory for DCM results
    output_dir = fullfile(spmeeg_path, 'dcm', 'M1', [input_type '_model']); 
    
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    prefix = 'bc_fl_ra_r_re_a_sm_ln_fh_spmeeg_';

    %% DCM 

    % Data filename
    %--------------------------------------------------------------------------
    filename = [prefix, p_name '_eeg.mat'];
    spmeeg_file = fullfile(spmeeg_path, filename);

    DCM.xY.Dfile = spmeeg_file;

    % Parameters and options used for setting up model
    %--------------------------------------------------------------------------
    DCM.options.analysis = 'ERP'; % analyze evoked responses
    DCM.options.model    = 'ERP'; % CMC model
    DCM.options.spatial  = 'ECD'; % spatial model
    DCM.options.trials   = [1, 2, 3, 4];     % index of ERPs within ERP/ERF file 
    DCM.options.Tdcm(1)  = 500;     % start of peri-stimulus time to be modelled
    DCM.options.Tdcm(2)  = dcm_end_time;   % end of peri-stimulus time to be modelled
    DCM.options.Nmodes   = 8;     % nr of modes for data selection 
    DCM.options.D        = 1;     % downsampling factor
    DCM.options.h        = 1;     % nr of DCT components %'1' for 'detrend' i.e. model the mean       
    DCM.options.onset    = 500 + 64;    % selection of onset in ms (prior mean); 64 is SPM default          
    DCM.options.dur      = 16;    % and dispersion (sd); 16ms is SPM default
    
    DCM.xY.modality = 'EEG'; % Specify modality 

    %--------------------------------------------------------------------------
    % Data and spatial model
    %--------------------------------------------------------------------------
    DCM  = spm_dcm_erp_data(DCM);
    
    %--------------------------------------------------------------------------
    % Location priors for dipoles
    %--------------------------------------------------------------------------
    DCM.Lpos  = indiv_mni_coords; 
    DCM.Sname = roi_labels;
    Nareas    = size(DCM.Lpos,2); %number of dipoles, used to create A/B/C matrices 
   
    %--------------------------------------------------------------------------
    % Spatial model
    %--------------------------------------------------------------------------
    DCM = spm_dcm_erp_dipfit(DCM);

    %--------------------------------------------------------------------------
    % Specify connectivity model
    %--------------------------------------------------------------------------
    
    %Setup forward and backward connection matrices
    DCM.A{1} = zeros(Nareas,Nareas); %forward matrix
    DCM.A{2} = zeros(Nareas,Nareas); %backward matrix
    DCM.A{3} = zeros(Nareas,Nareas); %lateral matrix
    
    DCM.C = zeros(Nareas); %setup empty C matrix

   % NOTE SPM USES (TARGET, SOURCE) FOR CONNECTIONS
    
    % Forward
    DCM.A{1}(3,1) = 1; % L_V1 -> L_FG
    DCM.A{1}(5,1) = 1; % L_V1 -> L_M1
    DCM.A{1}(7,1) = 1; % L_V1 -> L_PFC
    DCM.A{1}(4,2) = 1; % R_V1 -> R_FG
    DCM.A{1}(6,2) = 1; % R_V1 -> R_M1
    DCM.A{1}(8,2) = 1; % R_V1 -> R_PFC
    DCM.A{1}(5,3) = 1; % L_FG -> L_M1
    DCM.A{1}(7,3) = 1; % L_FG -> L_PFC
    DCM.A{1}(6,4) = 1; % R_FG -> R_M1
    DCM.A{1}(8,4) = 1; % R_FG -> R_PFC
    DCM.A{1}(7,5) = 1; % L_M1 -> L_PFC
    DCM.A{1}(8,6) = 1; % R_M1 -> R_PFC
    
    % Backward
    DCM.A{2}(1,3) = 1; % L_V1 <- L_FG
    DCM.A{2}(1,5) = 1; % L_V1 <- L_M1
    DCM.A{2}(1,7) = 1; % L_V1 <- L_PFC
    DCM.A{2}(2,4) = 1; % R_V1 <- R_FG
    DCM.A{2}(2,6) = 1; % R_V1 <- R_M1
    DCM.A{2}(2,8) = 1; % R_V1 <- R_PFC
    DCM.A{2}(3,5) = 1; % L_FG <- L_M1
    DCM.A{2}(3,7) = 1; % L_FG <- L_PFC
    DCM.A{2}(4,6) = 1; % R_FG <- R_M1
    DCM.A{2}(4,8) = 1; % R_FG <- R_PFC
    DCM.A{2}(5,7) = 1; % L_M1 <- L_PFC
    DCM.A{2}(6,8) = 1; % R_M1 <- R_PFC

    % Lateral
    DCM.A{3}(1,2) = 1; DCM.A{3}(2,1) = 1;  % V1
    DCM.A{3}(3,4) = 1; DCM.A{3}(4,3) = 1;  % FFA
    DCM.A{3}(5,6) = 1; DCM.A{3}(6,5) = 1;  % M1
    DCM.A{3}(7,8) = 1; DCM.A{3}(8,7) = 1;  % PFC
        
    % Modulation (B matrix)
    B_temp = zeros(Nareas);
    B_temp(3,1) = 1; B_temp(4,2) = 1;  % V1->FFA
    B_temp(5,1) = 1; B_temp(6,2) = 1;  % V1 -> M1
    B_temp(7,1) = 1; B_temp(8,2) = 1;  % V1 -> PFC
    B_temp(7,5) = 1; B_temp(8,6) = 1;  % M1 -> PFC
    B_temp(7,3) = 1; B_temp(8,4) = 1;  % FFA->PFC
    B_temp(5,3) = 1; B_temp(6,4) = 1;  % FFA->M1
    B_temp(1,3) = 1; B_temp(2,4) = 1;  % FFA->V1
    B_temp(3,7) = 1; B_temp(4,8) = 1;  % PFC->FFA
    B_temp(3,5) = 1; B_temp(4,6) = 1;  % FG<-M1
    B_temp(1,5) = 1; B_temp(2,6) = 1;  % V1 <- M1
    B_temp(1,7) = 1; B_temp(2,8) = 1;  % V1 <- PFC
    B_temp(5,7) = 1; B_temp(6,8) = 1;  % M1 <- PFC
    B_temp(1,2) = 1; B_temp(2,1) = 1;  % V1 lateral
    B_temp(3,4) = 1; B_temp(4,3) = 1;  % FFA lateral
    B_temp(5,6) = 1; B_temp(6,5) = 1;  % M1 lateral
    B_temp(7,8) = 1; B_temp(8,7) = 1;  % PFC lateral
    B_temp(1,1) = 1; B_temp(2,2) = 1;  % V1 self
    B_temp(3,3) = 1; B_temp(4,4) = 1;  % FFA self
    B_temp(5,5) = 1; B_temp(6,6) = 1;  % M1 self
    B_temp(7,7) = 1; B_temp(8,8) = 1;  % PFC self 

    DCM.B{1} = B_temp;  % Expectation
    % add more if you want to look at more effects 
    
    % expectation effect 
    DCM.xU.X = [-1;   % expected_neutral
                 1;   % unexpected_neutral
                -1;   % expected_fear
                 1];  % unexpected_fear

    DCM.xU.name = {'Expectation'};
    
    %DCM C matrix simply specified where the original input source is in the brain (V1 for this exp)
    DCM.C = [1; 1; 0; 0; 0; 0; 0; 0]; %Input sources (LV1 and RV1)
    
    % store meta-data
    DCM.input_type = input_type;
    DCM.M1 = 'M1 included';

    %--------------------------------------------------------------------------
    % Invert and save
    %--------------------------------------------------------------------------
    dcm_filename = ['DCM_M1_' input_type '_' p_name '.mat'];
    DCM.name = dcm_filename; 
    
    DCM = spm_dcm_erp(DCM); 

    disp(dcm_filename);
    save(fullfile(output_dir, dcm_filename), 'DCM');

    clear DCM Nareas dcm_filename spmeeg_file filename

end

