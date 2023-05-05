function fig = plot_connectivity(params,bidsID,connMeasure)

freqNames = fields(params.FreqBand)';    

% Load atlas
atlas = readtable(params.AtlasPath);
% Points by network
networks = {'Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
pos = cell(1,length(networks));
axisticks = ones(1,length(networks)+1);
axistickslabelpos = ones(1,length(networks));
for i=1:length(networks)
    pos{i}  = find(cellfun(@(x) contains(x,['_' networks{i} '_']), atlas.ROIName)); % all sources belonging to network{i}
    axisticks(i+1) = length(pos{i})+axisticks(i);
    axistickslabelpos(i) = axisticks(i)+(axisticks(i+1)-axisticks(i))/2; % set the network label in the middle
end
newpos = vertcat(pos{:});

fig = figure('Units','centimeters','Position',[0 0 10 8.5], 'visible', 'off');
tcl = tiledlayout(2,2,'TileSpacing','compact','Padding','none');

for iFreq=1:length(freqNames)
    try
        % Load connectivity matrices
        load(fullfile(params.ConnectivityPath,[bidsID '_' connMeasure '_' freqNames{iFreq} '.mat']),'connMatrix');
    catch
        error('Connectivity matrix could not be loaded');
    end
    
    % Reshape connectivity matrices to match the atlas networks
    connMatrix_r = connMatrix(newpos,newpos);
    n = size(connMatrix);
    
    % Plot
    ax = nexttile(tcl);
    im = imagesc(connMatrix_r);
    im.AlphaData = (triu(nan(n))+1); % Plot only half of the matrix as it is symetric
    set(ax,'XTick',axistickslabelpos,'XtickLabel',[],'YTick',axistickslabelpos,'YtickLabel',[],'TickLength',[0 0],'TickDir','out','box','off')
    if iFreq ==3
        set(ax,'XtickLabel',networks,'XtickLabelRotation',45,'YtickLabel',networks)
    end
    
    grid on
    set(ax,'GridColor','w','GridAlpha',0.5,'Layer','top');
    ax2 = copyobj(ax,ax.Parent);
    set(ax2,'Ytick',axisticks','Xtick',axisticks,'yticklabel',[],'xticklabel',[]);
    colorbar;
    title(freqNames(iFreq));
    
end
title(tcl,connMeasure);

end