% Parse params.json and set default values if not specified
%
% Cristina Gil, TUM, 01.02.2024
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

%% Checks to the params struct

% ====== PARALLELIZATION PARAMETERS =====
% NCores - OPTIONAL
if(~isfield(params,'NCores') || isempty(params.NCores))
    warning('NCores is not defined.\n')
elseif floor(params.NCores) ~= params.NCores
    error('NCores must be an integer.')
end


% ====== PATHS TO TOOLBOXES =======
% EEGLabPath - REQUIRED
if(~isfield(params,'EEGLabPath') || isempty(params.EEGLabPath) || (~isstring(params.EEGLabPath) && ~ischar(params.EEGLabPath)) || ~exist(fullfile(params.EEGLabPath,'eeglab.m'),'file'))
    error('Check EEGLabPath.');
end
% FieldTripPath - REQUIRED
if(~isfield(params,'FieldtripPath') || isempty(params.FieldtripPath) || (~isstring(params.FieldtripPath) && ~ischar(params.FieldtripPath)) || ~exist(fullfile(params.FieldtripPath,'ft_defaults.m'),'file'))
    error('Check FieldtripPath.');
end
% BrainConnectivityToolboxPath - REQUIRED
if(~isfield(params,'BrainConnectivityToolboxPath') || isempty(params.BrainConnectivityToolboxPath) || (~isstring(params.BrainConnectivityToolboxPath) && ~ischar(params.BrainConnectivityToolboxPath)))
    error('Check BrainConnectivityToolboxPath.');
end


% ====== DATA PARAMETERS ==========
% StudyName - OPTIONAL
if(~isfield(params,'StudyName') || isempty(params.StudyName))
    warning('StudyName is not defined.\n');
    params.StudyName = 'Study';
elseif ~isstring(params.StudyName) && ~ischar(params.StudyName)
    error('StudyName must be a string.');
end
% RawDataPath - REQUIRED
if(~isfield(params,'RawDataPath') || isempty(params.RawDataPath) || (~isstring(params.RawDataPath) && ~ischar(params.RawDataPath) )|| ~exist(params.RawDataPath,'dir'))
    error('Check RawDataPath.');
end
% PreprocessedDataPath - OPTIONAL
% If PreprocessedDataPath is not defined by the user, create a derivates folder with timestamp
if (~isfield(params,'PreprocessedDataPath') || isempty(params.PreprocessedDataPath))
    t = datestr(now,'yyyy_mm_dd');
    params.PreprocessedDataPath = fullfile(params.RawDataPath, ['derivatives_v' t]);
elseif (~isstring(params.PreprocessedDataPath) && ~ischar(params.PreprocessedDataPath))
    error('Check specefied PreprocessedDataPath.');
elseif ~exist(params.PreprocessedDataPath,'dir')
    mkdir(params.PreprocessedDataPath);
end
% Session - OPTIONAL
if (~isfield(params,'Session') || isempty(params.Session))
    params.Session = {};
elseif iscell(params.Session) 
elseif isstring(params.Session) || ischar(params.Session)
    params.Session = cellstr(params.Session);
else
    error('Session must be a string or array of strings.')
end
% Run - OPTIONAL
if (~isfield(params,'Run') || isempty(params.Run))
    params.Run = [];
elseif any(floor(params.Run) ~= params.Run)
    error('Run must be an integer or array of integers.')
end
% Task - OPTIONAL
if (~isfield(params,'Task') || isempty(params.Task))
    params.Task = [];
elseif ~ischar(params.Task)
    error('Task must be a string.')
end
% BidsChanloc - REQUIRED
if (~isfield(params,'BidsChanloc') || isempty(params.BidsChanloc) || ~any(strcmp(params.BidsChanloc,{'on','off'})))
    error('BidsChanloc must be "on" or "off".')
end
% NoseDir - OPTIONAL*
if (strcmp(params.BidsChanloc,'on') && (~isfield(params,'NoseDir') || isempty(params.NoseDir) || isempty(regexp(params.NoseDir,'[-+][XYZ]','once'))))
    error('If BidsChanloc is on, NoseDir has to be specified as +X, -X, +Y, -Y, +Z or -Z.')
end
% RefCoord - OPTIONAL*
if (strcmp(params.BidsChanloc,'on') && strcmp(params.AddRefChannel,'on') && ...
        (~isfield(params,'RefCoord') || isempty(params.RefCoord) || ...
        ~isfield(params.RefCoord,'X') || ~isfield(params.RefCoord,'Y') || ~isfield(params.RefCoord,'Z') || ...
        isempty(params.RefCoord.X) || isempty(params.RefCoord.Y) || isempty(params.RefCoord.Z) || ...
        ~isnumeric(params.RefCoord.X) || ~isnumeric(params.RefCoord.Y) || ~isnumeric(params.RefCoord.Z)))
    error('If BidsChanloc is on, and AddRefChannel is on, coordinates of the reference electrode need to be specified.')
end



% ====== PREPROCESSING PARAMETERS ======
% DownsamplingRate - OPTIONAL
if ~isfield(params,'DownsamplingRate')
    params.DownsamplingRate = [];
elseif floor(params.DownsamplingRate) ~= params.DownsamplingRate
    error('DownsamplingRate must be an integer.')
end
% FlatLineCriterion - REQUIRED
if ~isfield(params,'FlatLineCriterion') || isempty(params.FlatLineCriterion) || ...
        (~isnumeric(params.FlatLineCriterion) || (ischar(params.FlatLineCriterion) && ~strcmp(params.FlatLineCriterion,'off')))
    error("FlatLineCriterion must be specified as a number or set to 'off'");
end
% ChannelCriterion - REQUIRED
if ~isfield(params,'ChannelCriterion') || isempty(params.ChannelCriterion) || ...
        (~isnumeric(params.ChannelCriterion) || params.ChannelCriterion < 0 || params.ChannelCriterion > 1 || (ischar(params.ChannelCriterion) && ~strcmp(params.ChannelCriterion,'off')))
    error("ChannelCriterion must be specified as a number between 0 and 1 or set to 'off'");
end
% LineNoiseCriterion - REQUIRED
if ~isfield(params,'LineNoiseCriterion') || isempty(params.LineNoiseCriterion) || ...
        (~isnumeric(params.LineNoiseCriterion) || (ischar(params.LineNoiseCriterion) && ~strcmp(params.LineNoiseCriterion,'off')))
    error("LineNoiseCriterion must be specified as a number or set to 'off'");
end
% HighPass - REQUIRED
if ~isfield(params,'HighPass') || isempty(params.HighPass) || ...
        (any(~isnumeric(params.HighPass)) || length(params.HighPass) ~=2 || (ischar(params.HighPass) && ~strcmp(params.HighPass,'off')))
    error("HighPass must be specified as an array of two numbers or set to 'off'");
else
    params.HighPass = params.HighPass';
end
% AddRefChannel - REQUIRED
if (~isfield(params,'AddRefChannel') || isempty(params.AddRefChannel) || ~any(strcmp(params.AddRefChannel,{'on','off'})))
    error('AddRefChannel must be "on" or "off".')
end
% NICARepetitions - REQUIRED
if (~isfield(params,'NICARepetitions') || isempty(params.NICARepetitions)) || floor(params.NICARepetitions) ~= params.NICARepetitions || params.NICARepetitions <= 0
    error('NICARepetitions must be an integer higher than zero.')
end
% ICLabel - REQUIRED
if (~isfield(params,'ICLabel') || isempty(params.ICLabel)) || ~isnumeric(params.ICLabel) || any(size(params.ICLabel) ~= [7,2]) || ...
        any(any(params.ICLabel > 1)) || any(any(params.ICLabel < 0)) || any(params.ICLabel(:,2) < params.ICLabel(:,1))
    error('ICLabel must be a 7x2 array of max min probabilities for each component.')
end
% BurstCriterion - REQUIRED
if ~isfield(params,'BurstCriterion') || isempty(params.BurstCriterion) || ...
        (~isnumeric(params.BurstCriterion) || (ischar(params.BurstCriterion) && ~strcmp(params.BurstCriterion,'off')))
    error("BurstCriterion must be specified as a number or set to 'off'");
end
% WindowCriterion - REQUIRED
if ~isfield(params,'WindowCriterion') || isempty(params.WindowCriterion) || ...
        (~isnumeric(params.WindowCriterion) || params.WindowCriterion < 0 || params.WindowCriterion > 1 || (ischar(params.WindowCriterion) && ~strcmp(params.WindowCriterion,'off')))
    error("WindowCriterion must be specified as a number between 0 and 1 or set to 'off'");
end
% WindowCriterionTolerances - OPTIONAL
if ~isfield(params,'WindowCriterionTolerances') || isempty(params.WindowCriterionTolerances)
        warning("WindowCriterionTolerances was not specified. Using default")
        params.WindowCriterionTolerances = [];
else
    WinCrit = eval(params.WindowCriterionTolerances);
    if ~isnumeric(WinCrit) || length(WinCrit) ~=2
        error("WindowCriterionTolerances must be specified as an array of two numbers.");
    else
        params.WindowCriterionTolerances = WinCrit;
    end
end
% RejectBadTimeSegments - OPTIONAL
if ~isfield(params,'RejectBadTimeSegments') || isempty(params.RejectBadTimeSegments) 
    warning('RejectBadTimeSegments was not specified. Using default "on".');
    params.RejectBadTimeSegments = 'on';
elseif ~any(strcmp(params.RejectBadTimeSegments,{'on','off'}))
    error('RejectBadTimeSegments must be "on" or "off".')
end
% EpochLength - OPTIONAL
if ~isfield(params,'EpochLength') || isempty(params.EpochLength)
    params.EpochLength = [];
elseif ~isnumeric(params.EpochLength) || params.EpochLength <= 0
    error('EpochLength must be a number higher than zero.')
end
% EpochOverlap - OPTIONAL
if ~isfield(params,'EpochOverlap') || isempty(params.EpochOverlap)
    warning('EpochOverlap was not specified. Using default 0.');
    params.EpochOverlap = 0;
elseif ~isnumeric(params.EpochOverlap) || params.EpochOverlap < 0 || params.EpochOverlap > 1
    error('EpochLength must be a number between 0 and 1.')
end

% ====== EVENT-RELATED DATA PREPROCESSING PARAMETERS ======
% PreprocEventData - OPTIONAL
if ~isfield(params,'PreprocEventData') || isempty(params.BrainFeatExtr)
    params.PreprocEventData = false;    
end
if params.PreprocEventData == true
    disp('Event-related data will be preprocessed. Brain features will not be extracted.')
    params.BrainFeatExtr = false;
    if ~isempty(params.DownsamplingRate)
        params.DownsamplingRate = [];
        warning('Downsampling will not be carried out'); 
    end
    if ~strcmp(params.RejectBadTimeSegments,'off')
        params.RejectBadTimeSegments = 'off';
        warning('Bad time segments will not be rejected');
    end
    if ~isempty(params.EpochLength)
        params.EpochLength = [];
        warning('Data will not be segmented into epochs');
    end
end
% EventMarker - OPTIONAL
if params.PreprocEventData
    if (~isfield(params,'EventMarker') || isempty(params.EventMarker))
        error('EventMarker should be specified for preprocessing event-related data.');
    elseif ~ischar(params.EventMarker)
        error('EventMarker should be a string');
    end
end
% EventBounds - OPTIONAL
if params.PreprocEventData
    if (~isfield(params,'EventBounds') || isempty(params.EventBounds))
        error('EventBounds should be specified for preprocessing event-related data.')
    elseif ~isnumeric(params.EventBounds) || length(params.EventBounds) ~=2
        error("EventBounds must be specified as an array of two numbers.");
    end
end
    
    
% ====== FEATURE EXTRACTION PARAMETERS ======
% BrainFeatExtr - OPTIONAL
if ~isfield(params,'BrainFeatExtr') || isempty(params.BrainFeatExtr)
    params.BrainFeatExtr = true;
end
if params.BrainFeatExtr == false && params.PreprocEventData == false
    warning('Brain feature extraction will not be carried out because params.BrainFeatExtr was set to false by the user.')
elseif isempty(params.EpochLength) && params.PreprocEventData == false
    warning('Brain feature extraction will not be carried out because the parameter EpochLength was not specified.')
    params.BrainFeatExtr = false;
end
% FreqRes - OPTIONAL
if params.BrainFeatExtr
    if(~isfield(params,'FreqRes') || isempty(params.FreqRes))
        params.FreqRes = 1/EpochLength;
        warning('FreqRes was not specified. Using default 1/EpochLength: %0.3f',params.FreqRes);
    elseif ~isnumeric(params.FreqRes)
        error('FreqRes must be a number');
    elseif params.FreqRes < 1/params.EpochLength
        warning('FreqRes is smaller than 1/EpochLength. Zero padding will be used.');
        params.Pad = 1/params.FreqRes;
    end
else
    params.FreqRes = [];
end
% Pad - OPTIONAL
if params.BrainFeatExtr
    if ~isfield(params,'Pad')
        params.Pad = [];
    elseif ~isnumeric(params.Pad) || params.Pad <= params.EpochLength
        error('Pad must be a number higher than EpochLength.');
    elseif params.Pad ~= 1/params.FreqRes
        error('Pad is not the same as 1/FreqRes. Please specify one parameter or the other.')
    end
else
    params.Pad = [];
end
% FreqBand - OPTIONAL
if params.BrainFeatExtr
    if(~isfield(params,'FreqBand') || isempty(params.FreqBand))
        warning('FreqBand not specified. Using default frequency bands theta, alpha, beta and gamma');
        params.FreqBand.theta = [4 8-params.FreqRes];
        params.FreqBand.alpha = [8 13-params.FreqRes];
        params.FreqBand.beta = [13 30];
        params.FreqBand.gamma = [30+params.FreqRes 80];
    else
        bands = fields(params.FreqBand);
        for iBand = 1:length(bands)
            band = bands{iBand};
            if isempty(params.FreqBand.(band)) && strcmp(band,'theta')
                warning('Using default range for %s band: %d - %.3f,',band,4,8-params.FreqRes);
                params.FreqBand.theta = [4 8-params.FreqRes];
            elseif isempty(params.FreqBand.(band)) && strcmp(band,'alpha')
                warning('Using default range for %s band: %d - %.3f,',band,8,13-params.FreqRes);
                params.FreqBand.alpha = [8 13-params.FreqRes];
            elseif isempty(params.FreqBand.(band)) && strcmp(band,'beta')
                warning('Using default range for %s band: %d - %d,',band,13,30);
                params.FreqBand.beta = [13 30];
            elseif isempty(params.FreqBand.(band)) && strcmp(band,'gamma')
                warning('Using default range for %s band: %.3f - %d,',band,30+params.FreqRes,80);
                params.FreqBand.gamma = [30+params.FreqRes 80];
            elseif ~isempty(params.FreqBand.(band)) && isnumeric(params.FreqBand.(band)) ...
                    && length(params.FreqBand.(band)) == 2 && params.FreqBand.(band)(2) >= params.FreqBand.(band)(1)
                warning('Using user-specified range for %s band: %d - %d,',band,params.FreqBand.(band)(1),params.FreqBand.(band)(2));
            else
                error('Specify frequency bands as "band": [start end]')
            end
                
        end
    end
else
    params.FreqBand = [];
end
% Taper - OPTIONAL
if params.BrainFeatExtr
    if(~isfield(params,'Taper') || isempty(params.Taper))
        warning('Taper not specified. Using default "dpss".');
        params.Taper = 'dpss';
    end
else
    params.Taper = [];
end
% Tapsmofrq - OPTIONAL
if params.BrainFeatExtr
    if isfield(params,'Taper') && (~isfield(params,'Tapsmofrq') || isempty(params.Tapsmofrq))
        warning('Tapsmofrq not specified. Using default: 2');
        params.Tapsmofrq = 2;
    end
else
    params.Tapsmofrq = [];
end
% HeadModelPath - OPTIONAL
if params.BrainFeatExtr
    if ~isfield(params,'HeadModelPath') || isempty(params.HeadModelPath)
        warning('HeadModelPath not specified. Using default: standard_bem.mat');
        params.HeadModelPath = 'standard_bem.mat';
    end
else
    params.HeadModelPath = [];
end
% SurfaceModelPath - OPTIONAL
if params.BrainFeatExtr
    if ~isfield(params,'SurfaceModelPath') || isempty(params.SurfaceModelPath)
        warning('Surface<odelPath not specified. Using default: surface_white_both.mat');
        params.SurfaceModelPath = 'surface_white_both.mat';
    end
else
    params.SurfaceModelPath = [];
end
% AtlasPath - OPTIONAL
if params.BrainFeatExtr
    if ~isfield(params,'AtlasPath') || isempty(params.AtlasPath)
        warning('AtlasPath not specified. Using default: parcellations/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv');
        params.AtlasPath = 'parcellations/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';
    end
else
    params.AtlasPath = [];
end
% FreqResConnectivity - OPTIONAL
if params.BrainFeatExtr
    if ~isfield(params,'FreqResConnectivity') || isempty(params.FreqResConnectivity)
        warning('FreqResConnectivity not specified. Using default: 0.5');
        params.FreqResConnectivity = 0.5;
    elseif ~isnumeric(params.FreqResConnectivity) || params.FreqResConnectivity <=0
        error('FreqResConnectivity must be a number higher than zero.')
    end
else
    params.FreqResConnectivity = [];
end
% ConnMatrixThreshold - OPTIONAL
if params.BrainFeatExtr
    if ~isfield(params,'ConnMatrixThreshold') || isempty(params.ConnMatrixThreshold)
        warning('ConnMatrixThreshold not specified. Using default: 0.2');
        params.ConnMatrixThreshold = 0.2;
    elseif ~isnumeric(params.ConnMatrixThreshold) || params.ConnMatrixThreshold < 0 || params.ConnMatrixThreshold > 1
        error('ConnMatrixThreshold must be a number between 0 and 1.')
    end
else
    params.ConnMatrixThreshold = [];
end
%% Add toolboxes and pipeline functions to matlab path
if ~contains(path,params.EEGLabPath)
    run(fullfile(params.EEGLabPath,'eeglab.m'));
end
if ~contains(path,params.FieldtripPath)
    addpath(params.FieldtripPath);
    ft_defaults
end
addpath(params.BrainConnectivityToolboxPath);
addpath(fullfile('custom_functions'));
addpath(genpath(fullfile('external_functions')));

%% Create parpool object
if isempty(gcp('nocreate'))
    parpool(params.NCores);
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

if params.BrainFeatExtr
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
