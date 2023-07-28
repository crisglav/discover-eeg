% Generate plots for the figures 2 and 3 of the paper (peprocessing and EEG features)
%
% Cristina Gil Avila, 21.12.2022, Technical University of Munich

addpath ..
derivatives_path = '/rechenmagd3/Experiments/2021_preprocessing/datasets/LEMON-8min-closed-bids/derivatives_v2023_12_14/';
params = define_params(fullfile(derivatives_path,'params.json'));

% I moved the files since the preprocessing
params.PreprocessedDataPath = derivatives_path; 
params.PowerPath = fullfile(derivatives_path,'EEG_features','power');
params.SourcePath = fullfile(derivatives_path,'EEG_features','source');
params.ConnectivityPath = fullfile(derivatives_path,'EEG_features','connectivity');
params.GraphPath = fullfile(derivatives_path,'EEG_features','graph_measures');
params.AtlasPath = '../parcellations/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';
bidsID = 'sub-010002';
figures_path = '/rechenmagd3/Experiments/2021_preprocessing/figures';

%% FIGURE 2 - PREPROCESSING
bc_fig = plot_badchannels_singlestudy(params,bidsID);
saveas(bc_fig,fullfile(figures_path,[bidsID '_badchans.svg']));

ic_fig = plot_ICs_singlestudy(params,bidsID);
saveas(ic_fig,fullfile(figures_path,[bidsID '_ic.svg']));

bs_fig = plot_badtimesegments_singlestudy(params,bidsID);
saveas(bs_fig,fullfile(figures_path,[bidsID '_badsegments.svg']));

%% FIGURE 3 - FEATURE EXTRACTION
% Plot power and APF
[power_fig, topoplot_fig] = plot_power(params,bidsID);
saveas(power_fig,fullfile(figures_path,[bidsID '_power_apf.svg']));

% Plot source power
source_fig = plot_power_source(params,bidsID);
saveas(source_fig,fullfile(figures_path,[bidsID '_source.bmp']));

% Plot connectivity
dwpli_fig = plot_connectivity(params,bidsID,'dwpli');
saveas(dwpli_fig,fullfile(figures_path,[bidsID '_dwpli.svg']));
aec_fig = plot_connectivity(params,bidsID,'aec');
saveas(aec_fig,fullfile(figures_path,[bidsID '_aec.svg']));

% Plot graph measures
[dwpli_degree, dwpli_cc, dwpli_global]  = plot_graph_measures(params,bidsID,'dwpli');
saveas(dwpli_degree,fullfile(figures_path,[bidsID '_dwpli_degree.bmp']));
saveas(dwpli_cc,fullfile(figures_path,[bidsID '_dwpli_cc.bmp']));
saveas(dwpli_global,fullfile(figures_path,[bidsID '_dwpli_global.bmp']));

[aec_degree, aec_cc, aec_global]  = plot_graph_measures(params,bidsID,'aec');
saveas(aec_degree,fullfile(figures_path,[bidsID '_aec_degree.bmp']));
saveas(aec_cc,fullfile(figures_path,[bidsID '_aec_cc.bmp']));
saveas(aec_global,fullfile(figures_path,[bidsID '_aec_global.bmp']));

%%
f = figure('Units','centimeters','Position',[0 0 3 4]);
ax = axes;
c = colorbar(ax,'Location','southoutside');
caxis(dwpli_degree.Children(5).Children(1).Limits);
ax.Visible = 'off';
saveas(f,'colorbar_dwpli_degree.svg');

f = figure('Units','centimeters','Position',[0 0 3 4]);
ax = axes;
c = colorbar(ax,'Location','southoutside');
caxis(dwpli_cc.Children(5).Children(1).Limits);
ax.Visible = 'off';
saveas(f,'colorbar_cc_degree.svg');

f = figure('Units','centimeters','Position',[0 0 3 4]);
ax = axes;
c = colorbar(ax,'Location','southoutside');
caxis(aec_degree.Children(5).Children(1).Limits);
ax.Visible = 'off';
saveas(f,'colorbar_aec_degree.svg');

f = figure('Units','centimeters','Position',[0 0 3 4]);
ax = axes;
c = colorbar(ax,'Location','southoutside');
caxis(aec_cc.Children(5).Children(1).Limits);
ax.Visible = 'off';
saveas(f,'colorbar_aec_cc.svg');

