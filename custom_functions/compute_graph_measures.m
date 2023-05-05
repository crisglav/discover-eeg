function compute_graph_measures(params,bidsID,freqBand,connMeasure)
try
    load(fullfile(params.ConnectivityPath,[bidsID '_' connMeasure '_' freqBand '.mat']),'connMatrix');
catch
    error(['Connectivity matrix '  bidsID '_' connMeasure '_' freqBand ' could not be loaded']);
end

% Threshold the connectivity matrix depending on the desired amount of edges
sortedValues = sort(abs(connMatrix(:)),'descend');
sortedValues = sortedValues(~isnan(sortedValues));
threshold = sortedValues(floor(params.ConnMatrixThreshold*length(sortedValues)));

% Binarize connectivity matrix based on the threshold (adjacency matrix)
adjacency_matrix = abs(connMatrix) >= threshold;

%% Graph analysis measures from Brain Connectivity toolbox
% ---- Local measures ----- 
% Degree - Number of connexions of each node
degree = degrees_und(adjacency_matrix);
% Clustering coefficient - The percentage of existing triangles surrounding
% one node out of all possible triangles
cc = clustering_coef_bu(adjacency_matrix);

% ---- Global measures of segregation -----
% Global clustering coefficient
gcc = mean(cc);

% ---- Global measures of integration -----
% Characteristic path length
distance = distance_bin(adjacency_matrix);
[cpl, ~] = charpath(distance,0,0); % it does not include infinite distances in the calculation 
% Global efficiency - The average of the inverse shortest path between two points of the network
geff = efficiency_bin(adjacency_matrix); % it does include infinite distances in the calculation

% ----- Small-worldness -----
% Typically in small-world networks L >= Lrand but CC >> CCrand
randN = makerandCIJ_und(length(adjacency_matrix), floor(sum(sum(adjacency_matrix)))/2);
gcc_rand = mean(clustering_coef_bu(randN));
[cpl_rand, ~] = charpath(distance_bin(randN),0,0);
smallworldness = (gcc/gcc_rand) / (cpl/cpl_rand);

graph_measures.threshold = threshold;
graph_measures.degree = degree';
graph_measures.cc = cc;
graph_measures.gcc = gcc;
graph_measures.geff = geff;
graph_measures.smallworldness = smallworldness;
save(fullfile(params.GraphPath,[bidsID '_graph_' connMeasure '_' freqBand '.mat']),'graph_measures')
     
end
