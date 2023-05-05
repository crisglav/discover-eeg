function [f_degree, f_cc, f_global] = plot_graph_measures(params,bidsID,connMeasure)

freqNames = fields(params.FreqBand)';
graphMeas = {'gcc','geff','smallworldness'};
spider = nan(length(freqNames),length(graphMeas));

f_degree = figure('Units','centimeters','Position',[0 0 10 9], 'visible', 'off');
f_cc = figure('Units','centimeters','Position',[0 0 10 9], 'visible', 'off');
t_degree = tiledlayout(f_degree,2,2,'TileSpacing','none','Padding','compact');
t_cc = tiledlayout(f_cc,2,2,'TileSpacing','none','Padding','compact');

% Load surface structure
surf = ft_read_headshape('surface_white_both.mat');

% Load atlas positions
atlas = readtable(params.AtlasPath);
% Source model: centroid positons from Schaefer atlas
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = [atlas.R, atlas.A, atlas.S];
cfg.unit = 'mm';
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';
pos = sourcemodel_atlas.pos;
clear sourcemodel_atlas;

cmax_d = 0;
cmin_d = 0;
cmax_c = 0;
cmin_c = 0;
for iFreq=1:length(freqNames)
    
    % Load graph measures
    load(fullfile(params.GraphPath,[bidsID '_graph_' connMeasure '_' freqNames{iFreq} '.mat']),'graph_measures');
    
    % Plot degree at the 100 sources    
    % Interpolate colormap to correct range
    if  max(graph_measures.degree) > cmax_d
        cmax_d =  max(graph_measures.degree);
    end
    if  min(graph_measures.degree) < cmin_d
        cmin_d = min(graph_measures.degree);
    end
    index = fix((graph_measures.degree-cmin_d)/(cmax_d-cmin_d)*256)+1;
    rgb = squeeze(ind2rgb(index,parula(256)));
    
    ax_d(iFreq) = nexttile(t_degree);  
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.2);
    ft_plot_mesh(pos, 'vertexsize',12, 'vertexcolor',rgb);
    title(freqNames(iFreq));
    
    % Plot global clustering coefficient at the 100 sources
    if  max(graph_measures.degree) > cmax_c
        cmax_c =  max(graph_measures.cc);
    end
    if  min(graph_measures.degree) < cmin_c
        cmin_c = min(graph_measures.cc);
    end
    index = fix((graph_measures.cc-cmin_c)/(cmax_c-cmin_c)*256)+1;
    rgb = squeeze(ind2rgb(index,parula(256)));
    
    ax_c = nexttile(t_cc);  
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.2);
    ft_plot_mesh(pos, 'vertexsize',12, 'vertexcolor',rgb);
    title(freqNames(iFreq));

    % Set global graph measures in a matrix for plotting
    for iG = 1:length(graphMeas)
        spider(iFreq,iG) = graph_measures.(graphMeas{iG});
    end   
end

% title(t_degree,{connMeasure,'Degree'});
% title(t_cc,{connMeasure, 'Clust. coef.'});
% Set colormap and color limits for all subplots
set(ax_d, 'Colormap', parula, 'CLim', [cmin_d cmax_d])
set(ax_c, 'Colormap', parula, 'CLim', [cmin_c cmax_c])
% assign color bar to one tile
colorbar(ax_d(end),'eastoutside');
colorbar(ax_c(end),'eastoutside');

% Spider plot with global measures
f_global = figure('Units','centimeters','Position', [0 0 12 8], 'visible', 'off');
axeslim = [0, 0, min(spider(:,3)); 1, 1, max(spider(:,3))];
spider_plot(spider,...
    'AxesLabels', {'Global clustering coef.','Global efficiency','Smallworldness'}, ...
    'AxesLimits',axeslim, ...
    'AxesLabelsEdge','none',...
    'AxesFontSize',7,...
    'AxesAngular','off',...
    'MarkerSize',20);
legend(freqNames, 'Location', 'east');
ax_global = get(f_global,'CurrentAxes');
title(ax_global,[connMeasure ' - Global measures']);
legend(freqNames, 'Location', 'south');

end