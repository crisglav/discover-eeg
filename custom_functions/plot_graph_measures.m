function [f_degree, f_cc, f_global] = plot_graph_measures(params,bidsID,connMeasure)

freqNames = fields(params.freq_band)';
graphMeas = {'gcc','geff','smallworldness'};
spider = nan(length(freqNames),length(graphMeas));

f_degree = figure('Position',[412 412 1200 1200]);
f_cc = figure('Position',[412 412 1200 1200]);
t_degree = tiledlayout(f_degree,2,2);
t_cc = tiledlayout(f_cc,2,2);

% Load dummy source structure for plotting
load(fullfile(params.source_folder,[bidsID '_source_' freqNames{1} '.mat']),'source');

% Load surface structure
surf = ft_read_headshape('surface_white_both.mat');
for iFreq=1:length(freqNames)
    
    % Load graph measures
    load(fullfile(params.graph_folder,[bidsID '_graph_' connMeasure '_' freqNames{iFreq} '.mat']),'graph_measures');
        
    % Copy the pow structure and update with cc and degree
    source.avg.degree = graph_measures.degree;
    source.avg.cc = graph_measures.cc;

    % Interpolate local graph measures to the surface cortex
    cfg = [];
    cfg.method = 'nearest';
    cfg.parameter = {'degree','cc'};
    sourceInterp = ft_sourceinterpolate(cfg, source, surf);
       
    % Plot the interpolated data (Same as with ft_sourceplot but handling axes objects)
    degree = sourceInterp.degree;
    cc = sourceInterp.cc;
    
    nexttile(t_degree);
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv');
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', degree, 'clim', [min(degree) max(degree)],'colormap',jet(64));
    colorbar;
    camlight;  
    title(freqNames(iFreq));
    
    nexttile(t_cc);
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv');
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', cc, 'clim', [min(cc) max(cc)],'colormap',jet(64));
    colorbar;
    camlight;  
    title(freqNames(iFreq));

    % Set global graph measures in a matrix for plotting
    for iG = 1:length(graphMeas)
        spider(iFreq,iG) = graph_measures.(graphMeas{iG});
    end   
end

title(t_degree,{[connMeasure ' - Degree'],bidsID},'Interpreter','None');
title(t_cc,{[connMeasure ' - Clustering coef.'],bidsID},'Interpreter','None');

% Spider plot with global measures
f_global = figure('Position', [1988, 672, 780, 657]);
axeslim = [0, 0, min(spider(:,3)); 1, 1, max(spider(:,3))];
spider_plot(spider,...
    'AxesLabels', {'Global clustering coef.','Global efficiency','Smallworldness'}, ...
    'AxesLimits',axeslim, ...
    'AxesLabelsEdge','none');
ax_global = get(f_global,'CurrentAxes');
title(ax_global,{[connMeasure ' - Global measures'], bidsID},'Interpreter','None');
legend(freqNames, 'Location', 'south');

end