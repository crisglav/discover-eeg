function plot_electrodesandsources(params,bidsID)
% Plot volume conduction model, electrode positions and source positions

% Load one dataset to get electrode position
data = load_preprocessed_data(params,bidsID);

% Source model: centroid positons from Schaefer atlas
atlas400 = readtable(params.atlaspath);
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = [atlas400.R, atlas400.A, atlas400.S];
cfg.unit = 'mm';
cfg.headmodel = params.volpath;
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';

% Volume conduction model
load(params.volpath,'vol');

% Plotting
f = figure;
ft_plot_headmodel(vol,'facealpha',0.1,'facecolor',[0.1 0.1 0.1],'edgecolor',[1 1 1],'edgealpha',0.5);
hold on;
ft_plot_mesh(sourcemodel_atlas.pos, 'vertexsize',10, 'vertexcolor','b');
hold on;
ft_plot_sens(data.elec,'label','label','elec','true','elecshape','disc','elecsize',5,'facecolor','r');
view(90,0);
saveas(f,fullfile(params.figures_folder,'Electrodes and sources.svg'));
close(f)

end