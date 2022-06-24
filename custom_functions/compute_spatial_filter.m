function compute_spatial_filter(params,bidsID,freqBand)
% Load EEG data
data = load_preprocessed_data(params,bidsID);

%% Source model 
% Create a grid just for visualization (extract the positions outside the brain)
cfg = [];
cfg.method = 'basedonresolution';
cfg.resolution = 10;
cfg.unit = 'mm';
cfg.headmodel = params.volpath;
sourcemodel_grid = ft_prepare_sourcemodel(cfg);
sourcemodel_grid.coordsys = 'mni';
outside_pos = sourcemodel_grid.pos(~sourcemodel_grid.inside,:);

% Source model: centroid positons from Schaefer atlas
atlas400 = readtable(params.atlaspath);
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = cat(1,[atlas400.R, atlas400.A, atlas400.S],outside_pos);
cfg.sourcemodel.inside = [ones(size(atlas400,1),1); zeros(sum(~sourcemodel_grid.inside),1)];
cfg.unit = 'mm';
cfg.headmodel = params.volpath;
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';
% params.sourcemodel = sourcemodel_atlas;

%% Bandpass the data in the relevant frequency band
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = params.freq_band.(freqBand);
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
cfg.headmodel = params.volpath;
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

save(fullfile(params.source_folder,[bidsID '_source_' freqBand '.mat']),'source')

% %% Plot source power
% cfg = [];
% cfg.downsample = 2;
% cfg.parameter = 'pow';
% mri = ft_read_mri(params.mripath); % MRI template (for visualiztion)
% sourceInterp = ft_sourceinterpolate(cfg,source,mri);
% cfg = [];
% cfg.funparameter = 'pow';
% cfg.method = 'slice';
% ft_sourceplot(cfg,sourceInterp,mri);
% cfg = [];
% cfg.funparameter = 'pow';
% cfg.method = 'cloud';
% ft_sourceplot(cfg,source);

end