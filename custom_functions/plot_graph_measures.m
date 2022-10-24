function [f_degree, f_cc, f_global] = plot_graph_measures(params,bidsID,connMeasure)

freqNames = fields(params.freq_band)';
graphMeas = {'gcc','geff','smallworldness'};
spider = nan(length(freqNames),length(graphMeas));

f_degree = figure('Position',[412 412 1200 1200], 'visible', 'off');
f_cc = figure('Position',[412 412 1200 1200], 'visible', 'off');
t_degree = tiledlayout(f_degree,2,2);
t_cc = tiledlayout(f_cc,2,2);

% Load surface structure
surf = ft_read_headshape('surface_white_both.mat');

% Load atlas positions
atlas400 = readtable(params.atlaspath);
% Source model: centroid positons from Schaefer atlas
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = [atlas400.R, atlas400.A, atlas400.S];
cfg.unit = 'mm';
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';

for iFreq=1:length(freqNames)
    
    % Load graph measures
    load(fullfile(params.graph_folder,[bidsID '_graph_' connMeasure '_' freqNames{iFreq} '.mat']),'graph_measures');
    
    % Plot degree at the 400 sources
    % Interpolate colormap to correct range
    cmin = min(graph_measures.degree);
    cmax = max(graph_measures.degree);
    index = fix((graph_measures.degree-cmin)/(cmax-cmin)*256)+1;
    rgb = squeeze(ind2rgb(index,parula(256)));
    
    ax_d = nexttile(t_degree);  
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.1);
    ft_plot_mesh(sourcemodel_atlas.pos, 'vertexsize',20, 'vertexcolor',rgb);
    ax_d.CLim = [cmin cmax];
    colorbar(ax_d);
    title(freqNames(iFreq));
    
    % Plot global clustering coefficient at the 400 sources
    cmin = min(graph_measures.cc);
    cmax = max(graph_measures.cc);
    index = fix((graph_measures.cc-cmin)/(cmax-cmin)*256)+1;
    rgb = squeeze(ind2rgb(index,parula(256)));
    
    ax_c = nexttile(t_cc);  
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.1);
    ft_plot_mesh(sourcemodel_atlas.pos, 'vertexsize',20, 'vertexcolor',rgb);
    ax_c.CLim = [cmin cmax];
    colorbar(ax_c);
    title(freqNames(iFreq));

    % Set global graph measures in a matrix for plotting
    for iG = 1:length(graphMeas)
        spider(iFreq,iG) = graph_measures.(graphMeas{iG});
    end   
end

title(t_degree,{[connMeasure ' - Degree'],bidsID},'Interpreter','None');
title(t_cc,{[connMeasure ' - Clustering coef.'],bidsID},'Interpreter','None');

% Spider plot with global measures
f_global = figure('Position', [1988, 672, 780, 657], 'visible', 'off');
axeslim = [0, 0, min(spider(:,3)); 1, 1, max(spider(:,3))];
spider_plot(spider,...
    'AxesLabels', {'Global clustering coef.','Global efficiency','Smallworldness'}, ...
    'AxesLimits',axeslim, ...
    'AxesLabelsEdge','none');
ax_global = get(f_global,'CurrentAxes');
title(ax_global,{[connMeasure ' - Global measures'], bidsID},'Interpreter','None');
legend(freqNames, 'Location', 'south');

end