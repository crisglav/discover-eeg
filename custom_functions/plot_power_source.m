function fig = plot_power_source(params,bidsID)

freqNames = fields(params.FreqBand)';    

fig = figure('Units','centimeters','Position', [0 0 10 10],'Visible','off');
tcl = tiledlayout(2,2,'TileSpacing','compact','Padding','none');

% Load surface structure
surf = ft_read_headshape(params.SurfaceModelPath);

for iFreq=1:length(freqNames)
        
    % Load source file
    load(fullfile(params.SourcePath,[bidsID '_source_' freqNames{iFreq} '.mat']),'source');
    
    % Interpolate power to the surface cortex
    cfg = [];
    cfg.method = 'nearest';
    cfg.parameter = 'pow';
    sourceInterp = ft_sourceinterpolate(cfg, source, surf);
    
    % Plot the interpolated data (Same as with ft_sourceplot but handling axes objects)
    nexttile(tcl);  
    pow = sourceInterp.pow;
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv');
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', pow, 'clim', [min(pow) max(pow)],'colormap',parula);
    colorbar;
    camlight;
    title(freqNames(iFreq));
    
end

end