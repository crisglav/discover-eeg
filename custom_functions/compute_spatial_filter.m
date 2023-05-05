function compute_spatial_filter(params,bidsID,freqBand)
% Load EEG data
data = load_preprocessed_data(params,bidsID);

%% Source model 
% Source model: centroid positons from Schaefer atlas
atlas400 = readtable(params.AtlasPath);
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = [atlas400.R, atlas400.A, atlas400.S];
cfg.unit = 'mm';
cfg.headmodel = params.HeadModelPath;
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';

%% Bandpass the data in the relevant frequency band
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = params.FreqBand.(freqBand);
data = ft_preprocessing(cfg, data);

%% Compute the covariance matrix from the data
% First normalize time axis of the data (otherwise it cracks).
% Here we loose the temporal order of the epochs
temptime = data.time{1};
[data.time{:}] = deal(temptime);

% Compute the average covaraciance matrix from the sensor data
cfg = [];
cfg.covariance = 'yes';
cfg.keeptrials = 'no';
cfg.removemean = 'yes';
tlock = ft_timelockanalysis(cfg,data);

%%  Computation of the spatial filter
% Forward model (leadfield)
cfg = [];
cfg.sourcemodel = sourcemodel_atlas;
cfg.headmodel = params.HeadModelPath;
cfg.normalize = 'yes';
lf = ft_prepare_leadfield(cfg, data);

% Spatial filter
cfg = [];
cfg.method = 'lcmv';
cfg.keeptrials = 'yes';
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.lambda = '5%';
cfg.lcmv.fixedori = 'yes';
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.weightnorm = 'arraygain';
% cfg.lcmv.weightnorm = 'nai';
cfg.sourcemodel = lf;
source = ft_sourceanalysis(cfg, tlock);

save(fullfile(params.SourcePath,[bidsID '_source_' freqBand '.mat']),'source')

% Plot source power
% surf = ft_read_headshape('surface_white_both.mat');
% tmpcfg = [];
% tmpcfg.method = 'nearest';
% tmpcfg.parameter = 'pow';
% sourceInterp = ft_sourceinterpolate(tmpcfg, source, surf);
% 
% cfg = [];
% cfg.funparameter = 'pow';
% cfg.method = 'surface';
% cfg.funcolormap = 'jet';
% ft_sourceplot(cfg,sourceInterp,surf);

end