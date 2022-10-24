function fig = plot_power_source(params,bidsID)

freqNames = fields(params.freq_band)';    

fig = figure('Position',[412 412 1200 1200], 'visible', 'off');
tcl = tiledlayout(2,2);

% Load surface structure
surf = ft_read_headshape(params.surf);

for iFreq=1:length(freqNames)
    
    nexttile;
    
    % Load source file
    load(fullfile(params.source_folder,[bidsID '_source_' freqNames{iFreq} '.mat']),'source');

    % Interpolate power to the surface cortex
    cfg = [];
    cfg.method = 'nearest';
    cfg.parameter = 'pow';
    sourceInterp = ft_sourceinterpolate(cfg, source, surf);
    
    % Plot the interpolated data (Same as with ft_sourceplot but handling axes objects)
    pow = sourceInterp.pow;
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv');
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', pow, 'clim', [min(pow) max(pow)],'colormap',parula(64));
    colorbar;
    camlight;
   
    title(freqNames(iFreq));
    
end
title(tcl,{'Source power of ', bidsID},'Interpreter','None');

end