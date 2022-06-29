function f_global = plot_graph_measures(params,bidsID,connMeasure)
% Atlas networks - JUST FOR VISUALIZATION
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
axisticks = axisticks(2:end-1);

freqNames = fields(params.freq_band)';
graphMeas = params.graphMeas;
spider = nan(length(freqNames),length(graphMeas));

f_degree = figure('Position',[810   759   991   542]);
f_cc = figure('Position',[810   759   991   542]);
ax_degree = axes('Parent',f_degree);
ax_cc = axes('Parent',f_cc);

hold(ax_degree,'on')
hold(ax_cc,'on')
for iFreq=1:length(freqNames)
    
    % Load graph measures
    load(fullfile(params.graph_folder,[bidsID '_graph_' connMeasure '_' freqNames{iFreq} '.mat']),'graph_measures');
    
    % Local measures
    scatter(ax_degree,1:length(graph_measures.degree), graph_measures.degree,10,'filled');
    hold(ax_degree,'on')
    scatter(ax_cc,1:length(graph_measures.cc), graph_measures.cc,10,'filled');
    hold(ax_cc,'on');

    % Set global graph measures in a matrix for plotting
    for iG = 1:length(graphMeas)
        spider(iFreq,iG) = graph_measures.(graphMeas{iG});
    end   
end
% Draw vertical lines separating networks
for i=1:length(axisticks)
    xline(ax_degree,axisticks(i),'k','LineWidth',1)
    xline(ax_cc,axisticks(i),'k','LineWidth',1)
end
hold(ax_degree,'off');
hold(ax_cc,'off');
set(ax_degree,'XTick',axistickslabelpos,'XtickLabel',networks,'XtickLabelRotation',45,'TickLength',[0 0],'TickDir','out','box','off');
set(ax_cc,'XTick',axistickslabelpos,'XtickLabel',networks,'XtickLabelRotation',45,'TickLength',[0 0],'TickDir','out','box','off');
title(ax_cc,'Local clustering coefficient');
title(ax_degree,'Degree');
legend(ax_degree,freqNames);
legend(ax_cc,freqNames);
close(f_degree,f_cc) % NOT SAVED

axeslim = [0, 0, 0, 0, 0; 1, 1, max(spider(:,3)), 1, max(spider(:,5))];

f_global = figure;
spider_plot(spider,...
    'AxesLabels', graphMeas, ...
    'AxesLimits',axeslim);
title('Goblal graph measures')    
legend(freqNames, 'Location', 'south');

end