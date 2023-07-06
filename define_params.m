% Define preprocessing and brain features parameters
%
% Cristina Gil, TUM, 10.06.2022
function params = define_params(fname)

%% Read json file with all the parameters
if isempty(fname)
    fname = 'params_example.json';
end
fid = fopen(fname);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
params = jsondecode(str);

% If PreprocessedDataPath is not defined by the user, create a derivates
% folder with timestamp
if isempty(params.PreprocessedDataPath)
    t = datestr(now,'yyyy_mm_dd');
    params.PreprocessedDataPath = fullfile(params.RawDataPath, ['derivatives_v' t]);
end

% Convert '-Inf' to -Inf from params.json
WinCrit = params.WindowCriterionTolerances;
params.WindowCriterionTolerances = eval(params.WindowCriterionTolerances);

params.HighPass = params.HighPass';
%% Add toolboxes and pipeline functions to matlab path
run(fullfile(params.EEGLabPath,'eeglab.m'));
addpath(params.FieldtripPath);
ft_defaults
addpath(params.BrainConnectivityToolboxPath);
addpath(fullfile('custom_functions'));
addpath(genpath(fullfile('external_functions')));
addpath(pwd)

%% Create default parameters from params.json that are not specified by the user
% ===== Electrode positions =====
% Parameter to select electrode positions from the BIDS sidecar file (default 'off')
% If set to 'off' electrodes positions from a standard template in the MNI system are
% selected ('standard_1005.elc'). If you have recorded non-standard EEG channels or you have recorded
% the specific electrode positions of each participant and they are stored in a '*_electrodes.tsv' file,
% you can set this paramter to 'on' and the positions of the '*_electrodes.tsv' file will be used.
% Please if set to 'on', specify also the coordinate system of the recorded
% '*_electrodes.tsv' in params.NoseDir
if strcmp(params.BidsChanloc,'on')
    params.NoseDir = '+Y'; % 'RAS'
    % Note: Be careful with the orientation of the electrodes. If the coordinate system of the '*_electrodes.tsv'
    % is 'RAS' (see '*_coordsystem.json') you have to change the nose direction in EEGLAB to '+Y'.
else
    params.NoseDir = '+X'; % Default EEGLab coordinate system
end

% ===== Re-referencing =====
% Parameter to add back reference channel (default 'off')
% Set to 'on' if you want to add back the reference channel after rereferencing to the average reference. 
% In that case add also the coordinates of the reference channel in the same same coordinate system as electrodes.tsv or in the
% MNI system depending on th params.BidsChanloc
if strcmp(params.AddRefChannel,'on')
    % Coordinates of the reference electrode X, Y, Z in the same coordinate system as electrodes.tsv
    params.RefCoord.X = 0.3761; 
    params.RefCoord.Y = 27.39;
    params.RefCoord.Z = 88.668;
end

%  ==== Power spectrum ====
% Frequency resolution is by default 1/EpochLength, i.e. 0.5 Hz.
% It can be artificially increased with zero-padding. We select by default
% a zero-padding to have a frequency resolution of 0.1 Hz.
params.Padding = 10/params.EpochLength;

% ==== Frequency bands ==== 
% Frequency bands of interest (default values for theta, alpha, beta and gamma as in
% Pernet et al 2021)
if isempty(params.FreqBand.theta)
    params.FreqBand.theta = [4 8-params.FreqRes];
end
if isempty(params.FreqBand.alpha)
    params.FreqBand.alpha = [8 13-params.FreqRes];
end
if isempty(params.FreqBand.beta)
    params.FreqBand.beta = [13 30];
end
if isempty(params.FreqBand.gamma)
    params.FreqBand.gamma = [30+params.FreqRes 80];
end

%% Output folder structure
% Folder that will contain the preprocessing figures
params.FiguresPreprocessingPath = fullfile(params.PreprocessedDataPath,'preprocessing_visualization');
if ~exist(params.FiguresPreprocessingPath,'dir')
    mkdir(params.FiguresPreprocessingPath);
end
% Folder with reports
params.ReportsPath = fullfile(params.PreprocessedDataPath,'reports');
if ~exist(params.ReportsPath,'dir')
    mkdir(params.ReportsPath);
end

% Folder that will contain power files
params.PowerPath = fullfile(params.PreprocessedDataPath,'EEG_features','power');
if ~exist(params.PowerPath,'dir')
    mkdir(params.PowerPath);
end
% Folder that will contain source reconstruction files
params.SourcePath = fullfile(params.PreprocessedDataPath,'EEG_features','source');
if ~exist(params.SourcePath,'dir')
    mkdir(params.SourcePath);
end
% Folder that will contain connectivity matrices files
params.ConnectivityPath = fullfile(params.PreprocessedDataPath,'EEG_features','connectivity');
if ~exist(params.ConnectivityPath,'dir')
    mkdir(params.ConnectivityPath);
end
% Folder that will contain the graph-therory measures
params.GraphPath = fullfile(params.PreprocessedDataPath,'EEG_features','graph_measures');
if ~exist(params.GraphPath,'dir')
    mkdir(params.GraphPath);
end

%% Save params file to derivatives folder
% Convert -Inf to '-Inf'
paramsOut = params;
paramsOut.WindowCriterionTolerances = WinCrit;
jsonOutput = jsonencode(paramsOut);
fid = fopen(fullfile(paramsOut.PreprocessedDataPath,'params.json'),'w');
fprintf(fid,'%s',jsonOutput);
fclose(fid);

end
