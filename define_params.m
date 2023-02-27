% Define preprocessing and brain features parameters
%
% Cristina Gil, TUM, 10.06.2022
function params = define_params()
%% DATA AND TOOLBOXES PATHS (TO BE MODIFIED BY USER)

% The name of your study (use double quotes)
params.study = 'LONGITUDINAL';

% The path of your raw data in BIDS format
params.raw_data_path = '/rechenmagd4/Experiments/2023_1overf/data/Henrik';

% The path of the output of the pipeline (preprocessed data and extracted
% brain features). By default is stored in a created 'derivatives' folder.
t = datestr(now,'yyyy_mm_dd');
params.preprocessed_data_path = fullfile(params.raw_data_path, ['derivatives_v' t]);

% Parameters to select a specific task,run and and session.
params.session = {'baseline'}; % For all sessions: {} / for a specific session e.g. {'baseline'}
params.runs = [];
params.task =  'closed'; % For all tasks: [] / for a specific task e.g. 'closed'

% EEGLab path
params.eeglab_path = '/rechenmagd4/toolboxes_and_functions/eeglab';

% Fieldtrip path
params.fieldtrip_path = '/rechenmagd4/toolboxes_and_functions/fieldtrip';

% Brain Connectivity Toolbox path
params.bct_path = '/rechenmagd4/toolboxes_and_functions/2019_03_03_BCT';

% Add the toolboxes and pipeline functions to matlab path (does not need to be modified by the user)
run(fullfile(params.eeglab_path,'eeglab.m'));
addpath(params.fieldtrip_path);
ft_defaults
addpath(params.bct_path);
addpath(fullfile('custom_functions'));
addpath(genpath(fullfile('external_functions')));
addpath(pwd)
%% PREPROCESSING PARAMETERS
% ===== Electrode positions =====
% Parameter to select electrode positions from the BIDS sidecar file (default 'off')
% If set to 'off' electrodes positions from a standard template in the MNI system are
% selected ('standard_1005.elc'). If you have recorded non-standard EEG channels or you have recorded
% the specific electrode positions of each participant and they are stored in a '*_electrodes.tsv' file,
% you can set this paramter to 'on' and the positions of the '*_electrodes.tsv' file will be used.
% Please if set to 'on', specify also the coordinate system of the recorded
% '*_electrodes.tsv' in params.nosedir
params.bidschanloc = 'off'; 
if strcmp(params.bidschanloc,'on')
    params.nosedir = '+Y'; % 'RAS'
    % Note: Be careful with the orientation of the electrodes. If the coordinate system of the '*_electrodes.tsv'
    % is 'RAS' (see '*_coordsystem.json') you have to change the nose direction in EEGLAB to '+Y'.
else
    params.nosedir = '+X'; % Default EEGLab coordinate system
end
% Electrode template (Default 'standard_1005.elc')
% The default template has MNI coordinates and it is aligned with the
% standard BEM head model. If you have additional non-standard electrodes (e.g. LE and RE)
% you can define here your own electrode template (e.g. point to the '*_electrodes.tsv', but make sure that 
% it is properly aligned with the head model.
% params.elec_template = fullfile(params.fieldtrip_path,'template','electrode','standard_1005.elc');
% params.elec_template = '/rechenmagd4/Experiments/2021_preprocessing/datasets/CBP-mini/sub-CBPpa02/eeg/sub-CBPpa02_electrodes.tsv';

% ===== Downsampling ======
% In case you want to downsample the data, the new sampling date in Hz
params.sampling_rate = 250;

% ===== Bad channel rejection =====
% Parameters for the bad channel rejection (Default values as in clean_rawdata())
params.FlatlineCriterion = 5;
params.ChannelCriterion = 0.8;
params.LineNoiseCriterion = 4;
params.Highpass = [0.25 0.75];
params.FuseChanRej = 'off'; % If 'on' reject the union of bad channels detected in all the tasks.

% ===== Re-referencing =====
% Parameter to add back reference channel (default 'off')
% Set to 'on' if you want to add the reference channel as an active
% electrode after rereferencing to the average reference. In that case add
% also the coordinates of the reference channel in the same same coordinate system as electrodes.tsv or in the
% MNI system depending on th params.bidschanloc
params.addRefChannel = 'off';
if strcmp(params.addRefChannel,'on')
    % Coordinates of the reference electrode X, Y, Z in the same coordinate system as electrodes.tsv
    params.RefCoord.X = 0.3761; 
    params.RefCoord.Y = 27.39;
    params.RefCoord.Z = 88.668;
else
    params.RefCoord.X = [];
    params.RefCoord.Y = [];
    params.RefCoord.Z = [];
end

% ===== ICA =====
% Parameters for ICLabel (Default values as in Pernet et. al 2019,
% i.e. remove all ICs with probability of being muscle or eye > 80%)
params.IClabel = [nan nan; 0.8 1; 0.8 1; nan nan; nan nan; nan nan; nan nan]; 

% ===== Bad time segment removal =====
% Parameters for the bad time segments removal (Default values as in clean_rawdata()) 
params.BurstCriterion = 20;
params.WindowCriterion = 0.25;
params.WindowCriterionTolerances = [-Inf 7];

% ===== Segmentation into epochs =====
% Parameters to define the segmentation into epochs (Default 2s epochs with 50% overlap)
params.epoch_length = 2; % In seconds
params.epoch_overlap = 0.5; % In percentage

%%  FEATURE EXTRACTION PARAMETERS
%  ==== Power spectrum ====
% Frequency resolution is by default 1/epoch_length, i.e. 0.5 Hz.
% It can be artificially increased with zero-padding. We select by default
% a zero-padding to have a frequency resolution of 0.1 Hz.
params.padding = 10/params.epoch_length;
params.freq_res = 0.1;
% Type of taper for computing the power spectrum (default 'dpss')
params.taper = 'dpss';
% Frequency smoothing of the power spectrum (default 1)
params.tapsmofrq = 1;

% Frequency bands of interest (default values for theta, alpha, beta and gamma as in
% Pernet et al 2021
params.freq_band.theta = [4 8-params.freq_res];
params.freq_band.alpha = [8 13-params.freq_res]; 
params.freq_band.beta  = [13 30];
params.freq_band.gamma = [30+params.freq_res 80];             

% ==== Source reconstruction =====
% Head model (default 'standard_bem.mat')
params.volpath = fullfile(params.fieldtrip_path,'template','headmodel','standard_bem.mat');
% Atlas positions (default 100 source positions based on the Schaefer atlas)
params.atlaspath = fullfile(pwd,'parcellations','Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv');
% Brain surface model for plotting (default 'surface_white_both.mat')
params.surf = fullfile(params.fieldtrip_path,'template','anatomy','surface_white_both.mat');

% ===== Connectivity =====
% Resolution at which dwpli is computed (default 0.5),
% i.e. dwpli is computed each 0.5 Hz and averaged in the respective frequency band range
params.freq_res_connectivity = 0.5; % In Herz.

% ===== Graph measures =====
% Threshold at which the connectivity matrices are thersholded (default 0.2, i.e. 20% of top connections are kept)
% Several thresholds should be tested to check the robustness of the results.
params.connMatrix_threshold = 0.2; % In percentage

%% Generated folders (do not need to be modified by the user)
% Folder that will contain the preprocessing figures
params.figures_preprocessing_folder = fullfile(params.preprocessed_data_path,'preprocessing_visualization');
if ~exist(params.figures_preprocessing_folder,'dir')
    mkdir(params.figures_preprocessing_folder);
end
% Folder with reports
params.reports_folder = fullfile(params.preprocessed_data_path,'reports');
if ~exist(params.reports_folder,'dir')
    mkdir(params.reports_folder);
end

% Folder that will contain power files
params.power_folder = fullfile(params.preprocessed_data_path,'EEG_features','power');
if ~exist(params.power_folder,'dir')
    mkdir(params.power_folder);
end
% Folder that will contain source reconstruction files
params.source_folder = fullfile(params.preprocessed_data_path,'EEG_features','source');
if ~exist(params.source_folder,'dir')
    mkdir(params.source_folder);
end
% Folder that will contain connectivity matrices files
params.connectivity_folder = fullfile(params.preprocessed_data_path,'EEG_features','connectivity');
if ~exist(params.connectivity_folder,'dir')
    mkdir(params.connectivity_folder);
end
% Folder that will contain the graph-therory measures
params.graph_folder = fullfile(params.preprocessed_data_path,'EEG_features','graph_measures');
if ~exist(params.graph_folder,'dir')
    mkdir(params.graph_folder);
end

end
