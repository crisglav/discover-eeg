function graph_measures = compute_graph_measures(params,bidsID,freqBand,connMeasure)
try
    load(fullfile(params.connectivity_folder,[bidsID '_' connMeasure '_' freqBand '.mat']),'connMatrix');
catch
    error(['Connectivity matrix '  bidsID '_' connMeasure '_' freqBand ' could not be loaded']);
end

% Rearrange values of the connectivity matrix to match the atlas - JUST FOR VISUALIZATION
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
axisticks = axisticks(2:end);
newpos = vertcat(pos{:});
% Reshape connectivity matrices to match the atlas networks
c = connMatrix(newpos,newpos);

% Threshold the connectivity matrix depending on the desired amount of edges
sortedValues = sort(abs(connMatrix(:)),'descend');
sortedValues = sortedValues(~isnan(sortedValues));
threshold = sortedValues(floor(params.connMatrix_threshold*length(sortedValues)));

% Binarize connectivity matrix based on the threshold (adjacency matrix)
adjacency_matrix = abs(connMatrix) >= threshold;

%% Graph analysis measures from Brain Connectivity toolbox
% ---- Local measures ----- 
% Degree - Number of connexions of each node
degree = degrees_und(adjacency_matrix);

% Weighted degree - Number of connexions of each node wighted by the
% connexion value
% wdegree = nansum(adjacency_matrix.*connMatrix);

% Clustering coefficient - The percentage of existing triangles surrounding
% one node out of all possible triangles
clustering_coef = clustering_coef_bu(adjacency_matrix);

% ---- Global measures of segregation -----
% Global clustering coefficient
gcc = mean(clustering_coef);
% Transitivity
transitivity = transitivity_bu(adjacency_matrix);

% ---- Global measures of integration -----
% Characteristic path length
distance = distance_bin(adjacency_matrix);
[cpl, ~] = charpath(distance,0,0); % it does not include infinite distances in the calculation 
% Global efficiency - The average shortest path between two points of the network
geff = efficiency_bin(adjacency_matrix); % it does include infinite distances in the calculation

% ----- Small-worldness -----
% Typically in small-world networks L >= Lrand but CC >> CCrand
randN = makerandCIJ_und(length(adjacency_matrix), floor(sum(sum(adjacency_matrix)))/2);
gcc_rand = mean(clustering_coef_bu(randN));
[cpl_rand, ~] = charpath(distance_bin(randN),0,0);
smallworldness = (gcc/gcc_rand) / (cpl/cpl_rand);

graph_measures.threshold = threshold;
graph_measures.degree = degree;
graph_measures.cc = clustering_coef';
graph_measures.gcc = gcc;
graph_measures.transitivity = transitivity;
graph_measures.cpl = cpl;
graph_measures.geff = geff;
graph_measures.smallworldness = smallworldness;
save(fullfile(params.graph_folder,[bidsID '_graph_' connMeasure '_' freqBand '.mat']),'graph_measures')

% %% Graph analysis with BRAPH
% % g = GraphWU(connMatrix);
% g = GraphBU(connMatrix,'threshold',threshold);
% braph.degree = g.degree();
% % ---- Global measures of segregation -----
% [braph.globalcc, braph.localcc] = g.cluster();
% braph.transitivity = g.transitivity();
% % ---- Global measures of integration -----
% pth = g.pl();
% braph.cpl = g.measure(g.CPL_WSG);
% % ----- Small-worldness -----
% braph.smallworldness = smallworldness(g,1); % Takes a while


%% Plotting
% Thresholded connectivity matrix reorderd by networks
f1 = figure();
imagesc(c.*adjacency_matrix);
% imagesc(c,'AlphaData',adjacency_matrix+(~adjacency_matrix)*0.5);
ax = f1.CurrentAxes;
set(ax,'XTick',axistickslabelpos,'XtickLabel',networks,'XtickLabelRotation',45,'YTick',axistickslabelpos,'YtickLabel',networks,'TickLength',[0 0],'TickDir','out','box','off')
colorbar(ax);
grid on
grid minor
set(ax,'GridColor','w','GridAlpha',1,'LineWidth',1,'Layer','top','MinorGridColor','w','MinorGridLineStyle','-','MinorGridAlpha',0.2);
ax2 = copyobj(ax,ax.Parent);
set(ax2,'Ytick',axisticks','Xtick',axisticks,'yticklabel',[],'xticklabel',[],'Position',ax.Position);
ax2.XAxis.MinorTickValues = axistickslabelpos; % Subgrid separes left and right hemisferes
ax2.YAxis.MinorTickValues = axistickslabelpos;
% title(freqNames(iFreq));
    
    
figure()
scatter(1:400, graph_measures.degree,15,'filled');
for i=1:length(axisticks)
    hold on
    xline(axisticks(i),'k','LineWidth',1)
end
ax = gca;
set(ax,'XTick',axistickslabelpos,'XtickLabel',networks,'XtickLabelRotation',45,'TickLength',[0 0],'TickDir','out','box','off')
ylabel ('Number of connected nodes');
title('Degree');

% figure()
% scatter(1:400,wdegree,15,'filled')
% for i=1:length(axisticks)
%     hold on
%     xline(axisticks(i),'k','LineWidth',1)
% end
% ax = gca;
% set(ax,'XTick',axistickslabelpos,'XtickLabel',networks,'XtickLabelRotation',45,'TickLength',[0 0],'TickDir','out','box','off')
% ylabel ('Number of connected nodes weighted by its connectivity');
% title('Weighted degree')

f = figure();
scatter(1:400, graph_measures.clustering_coef,15,'filled');
for i=1:length(axisticks)
    hold on
    xline(axisticks(i),'k','LineWidth',1)
end
ax = f.CurrentAxes;
set(ax,'XTick',axistickslabelpos,'XtickLabel',networks,'XtickLabelRotation',45,'TickLength',[0 0],'TickDir','out','box','off')
ylabel ('Local clusterin coefficient'); % Fraction of nodes neighbours that are also neighbours of each other
title('Clustering coefficient');
% p = ax.InnerPosition;
% pos = [p(1)+p(1)*0.1,p(2), p(3)*0.5, p(4)*0.2];
% str = {'Global measures of segregation',sprintf('Global cc = %.2f', gcc),sprintf('Transitivity = %.2f', transitivity)};
% annotation('textbox',pos,'String',str,'BackgroundColor','w','FitBoxToText','on','verticalalignment','bottom');

end