function conn_fig = plot_connectivity(params,bidsID,connMeasure)

freqNames = fields(params.freq_band)';    

% Load atlas
atlas400 = readtable(params.atlaspath);
% Points by network
networks = {'Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
pos = cell(1,length(networks));
axisticks = ones(1,length(networks)+1);
axistickslabelpos = ones(1,length(networks));
for i=1:length(networks)
    pos{i}  = find(cellfun(@(x) contains(x,['_' networks{i} '_']), atlas400.ROIName)); % all sources belonging to network{i}
    axisticks(i+1) = length(pos{i})+axisticks(i);
    axistickslabelpos(i) = axisticks(i)+(axisticks(i+1)-axisticks(i))/2; % set the network label in the middle
end
newpos = vertcat(pos{:});

conn_fig = figure('Position',[412 412 1200 1200]);
tcl = tiledlayout(2,2);
for iFreq=1:length(freqNames)
    
    % Load connectivity matrices
    load(fullfile(params.connectivity_folder,[bidsID '_' connMeasure '_' freqNames{iFreq} '.mat']),'connMatrix');
    % Reshape connectivity matrices to match the atlas networks
    connMatrix_r = connMatrix(newpos,newpos);
%     clim = max(abs(max(max(connMatrix_r))),abs(min(min(connMatrix_r))));
    % Plot
    ax = nexttile;
    imagesc(connMatrix_r);
    set(ax,'XTick',axistickslabelpos,'XtickLabel',networks,'XtickLabelRotation',45,'YTick',axistickslabelpos,'YtickLabel',networks,'TickLength',[0 0],'TickDir','out','box','off')
    grid on
    grid minor
    set(ax,'GridColor','w','GridAlpha',1,'Layer','top','MinorGridColor','w','MinorGridLineStyle','-','MinorGridAlpha',0.2);
    ax2 = copyobj(ax,ax.Parent);
    set(ax2,'Ytick',axisticks','Xtick',axisticks,'yticklabel',[],'xticklabel',[]);
    ax2.XAxis.MinorTickValues = axistickslabelpos; % Subgrid separes left and right hemisferes
    ax2.YAxis.MinorTickValues = axistickslabelpos;
    colorbar;
    title(freqNames(iFreq));
    
end
title(tcl,[bidsID ' - Average ' connMeasure ],'Interpreter','None','Fontweight','bold');

% % Plot only one matrix at a freq band
% % Thresholded connectivity matrix reorderd by networks
% f1 = figure();
% imagesc(c.*adjacency_matrix);
% % imagesc(c,'AlphaData',adjacency_matrix+(~adjacency_matrix)*0.5);
% ax = f1.CurrentAxes;
% set(ax,'XTick',axistickslabelpos,'XtickLabel',networks,'XtickLabelRotation',45,'YTick',axistickslabelpos,'YtickLabel',networks,'TickLength',[0 0],'TickDir','out','box','off')
% colorbar(ax);
% grid on
% grid minor
% set(ax,'GridColor','w','GridAlpha',1,'LineWidth',1,'Layer','top','MinorGridColor','w','MinorGridLineStyle','-','MinorGridAlpha',0.2);
% ax2 = copyobj(ax,ax.Parent);
% set(ax2,'Ytick',axisticks','Xtick',axisticks,'yticklabel',[],'xticklabel',[],'Position',ax.Position);
% ax2.XAxis.MinorTickValues = axistickslabelpos; % Subgrid separes left and right hemisferes
% ax2.YAxis.MinorTickValues = axistickslabelpos;

end