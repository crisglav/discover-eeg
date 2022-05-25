% Define params
%
% 
function params = define_params()
%% Data paths (raw and preprocessed data folders in BIDS format)
params.study = 'TD-mini';
params.raw_data_path = '/rechenmagd4/Experiments/2021_preprocessing/TD-mini';
% params.raw_data_path = 'C:\Users\Mitarbeiter\eeg_datasets\rawBIDS';
params.preprocessed_data_path = fullfile(params.raw_data_path,'derivatives2');
params.task =  []; % 'restEC'; % [] % if you want all tasks


%% Toolboxes paths
% toolboxes_path = '/rechenmagd4/toolboxes_and_functions';
% EEGLab
params.eeglab_path = '/rechenmagd4/toolboxes_and_functions/eeglab';
% params.eeglab_path = 'C:\Users\Mitarbeiter\eeglab';
run(fullfile(params.eeglab_path,'eeglab.m'));
% % Clean_raw_data
% params.cleanrawdata_path = fullfile(toolboxes_path,'clean_rawdata-master');
% addpath(genpath(params.cleanrawdata_path));

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
addpath(fullfile('external_functions'));
%% PREPROCESSING PARAMETERS
% Electrode template (Use a standard one or yours)
params.elec_template = fullfile(params.fieldtrip_path,'template','electrode','standard_1005.elc');

% Bad channel rejection (Default parameters as defined in clean raw data)
params.FlatlineCriterion = 5;
params.ChannelCriterion = 0.8;
params.LineNoiseCriterion = 4;
params.Highpass = [0.25 0.75];

% ICLabel parameters
params.IClabel = [nan nan; 0.8 1; 0.8 1; nan nan; nan nan; nan nan; nan nan]; % Flag all ICs with probability of being muscle or eye > 80%

% Bad time windows removal (Default parameters as defined in clean raw data) 
params.BurstCriterion = 20;
params.WindowCriterion = 0.25;
params.WindowCriterionTolerances = [-Inf 7];

% Segmentation into epochs
params.epoch_length = 10; % In seconds
params.epoch_overlap = 0.5; % In percentage

%%  FEATURE EXTRACTION PARAMETERS
% Frequency resolution
params.freq_res = 1/params.epoch_length; % In Herz
params.freq_res_connectivity = 0.5; % In Herz. Increase the frequency resolution in computation of connectivity measures to speed up the computation time.

% Frequency bands of interest
params.freq_band.theta = [4 8-params.freq_res];  % Pernet et al 2021 Nat Neurosc (COBIDAS MEEG recomendations)
params.freq_band.alpha = [8 13-params.freq_res]; % Pernet et al 2021 Nat Neurosc (COBIDAS MEEG recomendations)
params.freq_band.beta = [13 30];                 % Pernet et al 2021 Nat Neurosc (COBIDAS MEEG recomendations)
params.freq_band.gamma = [60 100];               % Avoid line noise

% Peak frequency
params.taper = 'dpss';
params.tapsmofrq = 1;

% Source reconstruction
% Head model
params.volpath = fullfile(params.fieldtrip_path,'template','headmodel','standard_bem.mat');
% Atlas positions
params.atlaspath = fullfile('parcellations','Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv');

% Graph measures
params.connMatrix_threshold = 0.2;

%% Generated folders (don't need to change anything here)
% Folder that will contain the preprocessing figures
params.figures_preprocessing_folder = fullfile(params.preprocessed_data_path,'preprocessing_visualization');
if ~exist(params.figures_preprocessing_folder,'dir')
    mkdir(params.figures_preprocessing_folder);
end
% Folder that will contain the EEG_features figures
params.figures_folder = fullfile(params.preprocessed_data_path,'EEG_features','figures');
if ~exist(params.figures_folder,'dir')
    mkdir(params.figures_folder);
end

% Folder that will contain peak frequency files
params.pf_folder = fullfile(params.preprocessed_data_path,'EEG_features','peak_frequency');
if ~exist(params.pf_folder,'dir')
    mkdir(params.pf_folder);
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