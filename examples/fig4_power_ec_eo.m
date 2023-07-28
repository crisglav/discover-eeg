% Compare power spectrum between eyes open and eyes closed in the LEMON and
% TDBRAIN datasets
%
% Cristina Gil Avila, 19.07.2023, Technical University of Munich

clear all; close all;

% Add fieldtrip
addpath /rechenmagd4/toolboxes_and_functions/fieldtrip
ft_defaults
% Add raincloud plots
addpath(genpath('../external_functions/'));

figures_path = '/rechenmagd3/Experiments/2021_preprocessing/figures';
results_path = '/rechenmagd3/Experiments/2021_preprocessing/results';

conditions = {'closed','open'};
colors = [0 0.447 0.741; 0.455 0.674 0.188];
power_fig = figure('Units','centimeters','Position', [0 0 18 7],'Visible','on');
tlc = tiledlayout(1,2);
nexttile
%% LEMON dataset - list recordings of participans with eyes open and closed

% List recordings of participans with eyes open and closed
for c=1:length(conditions)
    condition = conditions{c};
    switch condition
        case 'closed'
            study_path = '/rechenmagd3/Experiments/2021_preprocessing/datasets/LEMON-8min-closed-bids/derivatives_v2023_05_07';
        case 'open'
            study_path = '/rechenmagd3/Experiments/2021_preprocessing/datasets/LEMON-8min-open-bids/derivatives_v2023_05_05';
    end
    % Load power spectrum from all recordings
    power_files_lemon.(condition) = dir(fullfile(study_path, 'EEG_features','power','*_power.mat'));
end
cl = {power_files_lemon.closed.name};
op = {power_files_lemon.open.name};
noop = find(~ismember(cl,op));
nocl = find(~ismember(op,cl));
% Note that there is one participant that has an eyes open recording only
% ('sub-010126_power.mat'). This is because the marker in the original data
% of the eyes closed condition was S208 instead of S210. We will discard
% this participant from the analysis.

% Load power spectra
for c=1:length(conditions)
    
    condition = conditions{c};  
    nSubj = length(power_files_lemon.(condition));
    aux = cell(1,nSubj);
    % Load power spectrum
    for i=1:nSubj
        if strcmp(condition,'open') && i==nocl % Skip 'sub-010126_power.mat'
            continue;
        end
        load(fullfile(power_files_lemon.(condition)(i).folder,power_files_lemon.(condition)(i).name));
        aux{i} = power;
    end
          
    aux = aux(~cellfun('isempty',aux));
    participants_lemon.(condition) = aux;
    clear aux;
end

l=cell(1,length(conditions));

for c=1:length(conditions)
    condition = conditions{c};
    % Grand average across participants
    cfg = [];
    cfg.channel = 'all';
    cfg.parameter = 'powspctrm';
    cfg.keepindividual = 'yes';
    grandavg_lemon.(condition) = ft_freqgrandaverage(cfg,participants_lemon.(condition){:});
       
%     % Plot Pz - To compare with van Dijk et al. 2022 Figure 4
%     cfg = [];
%     cfg.channel = 'Pz';
%     pz = ft_selectdata(cfg, grandavg_lemon.(condition));
%     
%     l{c} = stdshade(log(squeeze(pz.powspctrm)),0.1,colors(c,:),pz.freq);
    
    % Plot channel average
    cfg = [];
    cfg.channel = 'all';
    cfg.avgoverchan = 'yes';
    avg = ft_selectdata(cfg, grandavg_lemon.(condition));
    
    l{c} = stdshade(log(squeeze(avg.powspctrm)),0.1,colors(c,:),avg.freq);
    hold on
end

title('LEMON Power spectrum (chann avg)');
ylabel('logPower');
xlabel('Frequency (Hz)');
ylim([-8 2]);
box off

% Statistical analysis: dependent samples cluster based permutation test across frequencies
nSubj = length(cl);

cfg = [];
% cfg.channel = 'Pz';
cfg.avgoverchan = 'yes';
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.parameter = 'powspctrm';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.alpha = 0.025;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.numrandomization = 500;
design = [1:nSubj, 1:nSubj;ones(1,nSubj), 2*ones(1,nSubj)];
cfg.design = design;
cfg.uver = 1; % Dependent variable: participant number
cfg.ivar = 2; % Independent variable: Condition. EC = 1, EO = 2
cfg.avgoverfreq = 'no';

stat_lemon = ft_freqstatistics(cfg, participants_lemon.closed{:}, participants_lemon.open{:});
save(fullfile(results_path,'stats_power_lemon.mat'),'stat_lemon');

% Get relevant values
pos_cluster_pvals = [stat_lemon.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos = ismember(stat_lemon.posclusterslabelmat, pos_clust);

neg_cluster_pvals = [stat_lemon.negclusters(:).prob];
neg_clust = find(neg_cluster_pvals < 0.025);
neg = ismember(stat_lemon.negclusterslabelmat, neg_clust);

significant_freqs = grandavg_lemon.closed.freq(pos);
significant_freqs_neg = grandavg_lemon.closed.freq(neg);

plot(significant_freqs,-8*ones(1,length(significant_freqs)),'k','LineWidth',3);
plot(significant_freqs_neg,-8*ones(1,length(significant_freqs_neg)),'k','LineWidth',3);
legend([l{:}],conditions);

%% TDBRAIN dataset
study_path = '/rechenmagd3/Experiments/2021_preprocessing/datasets/vanDijk/derivatives_v2023_07_14';
files_closed = dir(fullfile(study_path, 'EEG_features','power','*task-restEC_power.mat'));
files_open = dir(fullfile(study_path, 'EEG_features','power','*task-restEO_power.mat'));

% Select only one recording per participants (some participants have more
% than one session)
ids_ec = cellfun(@(x) x(1:12),{files_closed.name},'UniformOutput',false);
[~,ic] = unique(ids_ec);
power_files_vd.closed = files_closed(ic);

ids_eo = cellfun(@(x) x(1:12),{files_open.name},'UniformOutput',false);
[~,io] = unique(ids_eo);
power_files_vd.open = files_open(io);

nSubj_vd = length(ic);

% Load power spectra
for c=1:length(conditions)
    
    condition = conditions{c};  
    aux = cell(1,nSubj_vd);
    % Load power spectrum
    for i=1:nSubj_vd
        load(fullfile(power_files_vd.(condition)(i).folder,power_files_vd.(condition)(i).name));
        aux{i} = power;
    end
    aux = aux(~cellfun('isempty',aux));
    participants_vd.(condition) = aux;
    clear aux;
end

% Perform grand average, and plot average power spectra in eyes open and closed
nexttile
l=cell(1,length(conditions));
for c=1:length(conditions)
    condition = conditions{c};  
    % Grand average across participants
    cfg = [];
    cfg.channel = 'all';
    cfg.parameter = 'powspctrm';
    cfg.keepindividual = 'yes';
    grandavg_vd.(condition) = ft_freqgrandaverage(cfg,participants_vd.(condition){:});
       
%     % Plot Pz  - To compare with van Dijk et al. 2022 Figure 4
%     cfg = [];
%     cfg.channel = 'Pz';
%     pz = ft_selectdata(cfg, grandavg_vd.(condition));
%     l{c} = stdshade(log(squeeze(pz.powspctrm)),0.1,colors(c,:),pz.freq);
    
    % Plot channel average
    cfg = [];
    cfg.channel = 'all';
    cfg.avgoverchan = 'yes';
    avg = ft_selectdata(cfg, grandavg_vd.(condition));
    l{c} = stdshade(log(squeeze(avg.powspctrm)),0.1,colors(c,:),avg.freq);

    hold on
end

title('TDBRAIN Power spectrum (chan avg)');
ylabel('logPower');
xlabel('Frequency (Hz)');
ylim([-8 2]);
box off

% Statistical analysis: dependent samples cluster based permutation test across frequencies
cfg = [];
% cfg.channel = 'Pz';
cfg.avgoverchan = 'yes';
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.parameter = 'powspctrm';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.alpha = 0.025;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.numrandomization = 500;
design = [1:nSubj_vd, 1:nSubj_vd;ones(1,nSubj_vd), 2*ones(1,nSubj_vd)];
cfg.design = design;
cfg.uver = 1; % Dependent variable: participant number
cfg.ivar = 2; % Independent variable: Condition. EC = 1, EO = 2
cfg.avgoverfreq = 'no';

stat_vd = ft_freqstatistics(cfg, participants_vd.closed{:}, participants_vd.open{:});
save(fullfile(results_path,'stats_power_vanDijk.mat'),'stat_vd');
% Get relevant values
pos_cluster_pvals = [stat_vd.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos = ismember(stat_vd.posclusterslabelmat, pos_clust);

neg_cluster_pvals = [stat_vd.negclusters(:).prob];
neg_clust = find(neg_cluster_pvals < 0.025);
neg = ismember(stat_vd.negclusterslabelmat, neg_clust);

significant_freqs = grandavg_vd.closed.freq(pos);
significant_freqs_neg = grandavg_vd.closed.freq(neg);

% Plot a black line indicating where power is significantly different
plot(significant_freqs,-8*ones(1,length(significant_freqs)),'k','LineWidth',3);
plot(significant_freqs_neg,-8*ones(1,length(significant_freqs_neg)),'k','LineWidth',3);
legend([l{:}],conditions);

saveas(power_fig,fullfile(figures_path,'power_lemon_vd_avg.svg'));