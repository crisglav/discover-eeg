% Define preprocessing and brain features parameters
%
% Cristina Gil, TUM, 10.06.2022
function params = define_params()
%% Data paths (raw and preprocessed data folders in BIDS format)
params.study = 'CBP-mini';
params.raw_data_path = '/rechenmagd4/Experiments/2021_preprocessing/datasets/CBP-mini'; 
% params.raw_data_path = '/rechenmagd3/CGX_MS/EEG/rawBIDS';
% params.raw_data_path = 'C:\Users\Mitarbeiter\eeg_datasets\rawBIDS';
params.preprocessed_data_path = fullfile(params.raw_data_path, 'derivatives'); 
% params.preprocessed_data_path = '/rechenmagd3/CGX_MS/EEG/preprocessed_Pernet2019_Interpolationv4';
params.task =  []'; % 'restEC'; % [] % if you want all tasks


%% Toolboxes paths
% EEGLab
params.eeglab_path = '/rechenmagd4/toolboxes_and_functions/eeglab';
% params.eeglab_path = 'C:\Users\Mitarbeiter\eeglab';
run(fullfile(params.eeglab_path,'eeglab.m'));

% Fieldtrip
params.fieldtrip_path = '/rechenmagd4/toolboxes_and_functions/fieldtrip'; % fresh from github
% params.fieldtrip_path = 'C:\Users\Mitarbeiter\fieldtrip'; % fresh from github
addpath(params.fieldtrip_path);
ft_defaults

% Brain Connectivity Toolbox
params.bct_path = '/rechenmagd4/toolboxes_and_functions/2019_03_03_BCT';
% params.bct_path = 'C:\Users\Mitarbeiter\bct\2019_03_03_BCT';
addpath(params.bct_path);

% BRAPH
params.braph = '/rechenmagd4/toolboxes_and_functions/BRAPH/BRAPH 1.0.0/graph';
addpath(params.braph);

% Custom functions
addpath(fullfile('custom_functions'));
% External functions (mainly from matlab central)
% addpath(fullfile('external_functions'));
%% PREPROCESSING PARAMETERS
params.bidschanloc = 'on'; % Default 'off', only set to on if you have nonstandard EEG channels or you recorded the electrodes positions
params.nosedir = '+Y'; % 'RAS' % Needs to be specified in case bidschanloc = 'on'
% Note: if you want to use the electrode location of the electrodes.tsv
% file you can change 'bidschanloc to 'on'. Be careful with the orientation
% of the electrodes. If the coordinate system of the electrodes.tsv is 'RAS'
% (See *_coordsystem.json) you have to change the nose direction in EEGLAB to '+Y'
% Note2: bidschanloc 'on' does not work if one specific task is selected. 

% Add back reference channel
params.addRefChannel = 'yes';
% Coordinates of the reference electrode X, Y, Z in the same coordinate system as electrodes.tsv
params.RefCoord.X = 0.3761; 
params.RefCoord.Y = 27.39;
params.RefCoord.Z = 88.668;
if strcmp(params.addRefChannel,'no') || strcmp(params.addRefChannel,'off')
    params.RefCoord.X = [];
    params.RefCoord.Y = [];
    params.RefCoord.Z = [];

end

% Bad channel rejection (Default parameters as defined in clean_rawdata())
params.FlatlineCriterion = 5;
params.ChannelCriterion = 0.8;
params.LineNoiseCriterion = 4;
params.Highpass = [0.25 0.75];
params.FuseChanRej = 'off'; % If 'on' reject the union of bad channels detected in all the tasks.

% ICLabel parameters (Default: Flag all ICs with probability of being muscle or eye > 80%)
params.IClabel = [nan nan; 0.8 1; 0.8 1; nan nan; nan nan; nan nan; nan nan]; 

% Bad time windows removal (Default parameters as defined in clean_rawdata()) 
params.BurstCriterion = 20;
params.WindowCriterion = 0.25;
params.WindowCriterionTolerances = [-Inf 7];

% Segmentation into epochs
params.epoch_length = 10; % In seconds
params.epoch_overlap = 0.5; % In percentage

%%  FEATURE EXTRACTION PARAMETERS
% Electrode template (Use a standard one or yours)
% This template has MNI coordinates and it is aligned with the standard BEM
% model. If you have additional non-standard electrodes (e.g. LE and RE)
% you can define here your electrode template, but make sure that it is
% properly aligned with the head model.
params.elec_template = fullfile(params.fieldtrip_path,'template','electrode','standard_1005.elc');
% params.elec_template = '/rechenmagd4/Experiments/2021_preprocessing/datasets/CBP-mini/sub-CBPpa02/eeg/sub-CBPpa02_electrodes.tsv';

% Frequency resolution
params.freq_res = 1/params.epoch_length; % In Herz
params.freq_res_connectivity = 0.5; % In Herz. Increase the frequency resolution in computation of connectivity measures to speed up the computation time.

% Frequency bands of interest
params.freq_band.theta = [4 8-params.freq_res];  % Pernet et al 2021 Nat Neurosc (COBIDAS MEEG recomendations)
params.freq_band.alpha = [8 13-params.freq_res]; % Pernet et al 2021 Nat Neurosc (COBIDAS MEEG recomendations)
params.freq_band.beta = [13 30];                 % Pernet et al 2021 Nat Neurosc (COBIDAS MEEG recomendations)
params.freq_band.gamma = [60 100];               % Avoid line noise at 50 Hz

% Peak frequency
% For the computation of the power spectrum a multitaper method is used
params.taper = 'dpss';
params.tapsmofrq = 1;

% Source reconstruction
% Head model is standard BEM
params.volpath = fullfile(params.fieldtrip_path,'template','headmodel','standard_bem.mat');
% Atlas positions
params.atlaspath = fullfile('parcellations','Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv');
params.mripath = fullfile(params.fieldtrip_path,'template','anatomy','single_subj_T1.nii');

% Graph measures
params.connMatrix_threshold = 0.2;

%% Generated folders (don't need to change anything here)
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
% Folder that will contain the EEG_features figures
params.figures_folder = fullfile(params.preprocessed_data_path,'EEG_features','figures');
if ~exist(params.figures_folder,'dir')
    mkdir(params.figures_folder);
end

% % Folder that will contain peak frequency files
% params.pf_folder = fullfile(params.preprocessed_data_path,'EEG_features','peak_frequency');
% if ~exist(params.pf_folder,'dir')
%     mkdir(params.pf_folder);
% end
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