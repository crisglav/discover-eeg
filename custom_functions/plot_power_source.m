function fig = plot_power_source(params,bidsID)

freqNames = fields(params.freq_band)';    

fig = figure('Position',[412 412 1200 1200], 'visible', 'off');
tcl = tiledlayout(2,2);
tcl.TileSpacing = 'compact';
tcl.Padding = 'compact';
% Load surface structure
surf = ft_read_headshape(params.surf);

for iFreq=1:length(freqNames)
        
    % Load source file
    load(fullfile(params.source_folder,[bidsID '_source_' freqNames{iFreq} '.mat']),'source');
    
    % Plot power estimates at 400 sources
%     % Interpolate colormap to correct range
%     cmin = min(source.avg.pow);
%     cmax = max(source.avg.pow);
%     index = fix((source.avg.pow-cmin)/(cmax-cmin)*256)+1;
%     rgb = squeeze(ind2rgb(index,parula(256)));
%     
%     ax = nexttile(tcl);  
%     ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.1);
%     ft_plot_mesh(source.pos, 'vertexsize',20, 'vertexcolor',rgb);
%     ax.CLim = [cmin cmax];
%     c = colorbar(ax);
%     c.Label.String = 'Power (uV^2/Hz)';
%     title(freqNames(iFreq));
    
    
    % Interpolate power to the surface cortex
    cfg = [];
    cfg.method = 'nearest';
    cfg.parameter = 'pow';
    sourceInterp = ft_sourceinterpolate(cfg, source, surf);
    
    % Plot the interpolated data (Same as with ft_sourceplot but handling axes objects)
    nexttile(tcl);  
    pow = sourceInterp.pow;
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv');
    ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', pow, 'clim', [min(pow) max(pow)],'colormap',parula(64));
    colorbar;
    camlight;
    title(freqNames(iFreq));
    
end
title(tcl,{'Source power of ', bidsID},'Interpreter','None');

end