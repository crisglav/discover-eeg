function plot_atlasregions(params)
% Source model: centroid positons from Schaefer atlas
atlas = readtable(params.AtlasPath);
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = [atlas.R, atlas.A, atlas.S];
cfg.unit = 'mm';
cfg.headmodel = params.HeadModelPath;
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';

% Volume conduction model
load(params.HeadModelPath,'vol');

% Points by network
networks = {'Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
c = lines(length(networks));
f = figure('Position',[651 613 972 684]);
ft_plot_mesh(vol.bnd(3),'facealpha',0.1,'facecolor',[0.1 0.1 0.1],'edgecolor',[1 1 1],'edgealpha',0.5); % brain
for i=1:length(networks)
    hold on
    idx  = find(cellfun(@(x) contains(x,['_' networks{i} '_']), atlas.ROIName));
    ft_plot_mesh(sourcemodel_atlas.pos(idx,:), 'vertexsize',15, 'vertexcolor',c(i,:));
    h(i) = scatter(NaN,NaN,10,c(i,:),'filled'); % dummy for legend
end
legend(h,networks,'color','none');
view(90,0);
title ('Schaefer atlas 100 sources')
saveas(f,fullfile(params.preprocessed_data_path,'EEG_features','Schaefer atlas.svg'));
close(f);
end