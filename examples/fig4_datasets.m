%% Datasets demographics and preprocessing.
%
% Code to generate first part of figure 4 of the paper.
%
% Cristina Gil, 28.07.2023, Technical University of Munich

% Add toolboxes to the path
clear all, close all; 
addpath ../external_functions/raincloudplots
addpath ../custom_functions
addpath ..

figures_path = '/rechenmagd3/Experiments/2021_preprocessing/figures';
results_path = '/rechenmagd3/Experiments/2021_preprocessing/results';
%% Datasets
% Note that there is one participant in the LEMON study that has an eyes open recording only
% ('sub-010126_power.mat'). This is because the marker in the original data
% of the eyes closed condition was S208 instead of S210. We will discard
% this participant from the analysis.

% LEMON eyes closed dataset
derivatives_path = '/rechenmagd3/Experiments/2021_preprocessing/datasets/LEMON-8min-closed-bids/derivatives_v2023_05_07';
params_lemon = define_params(fullfile(derivatives_path,'params.json'));
params_lemon.PreprocessedDataPath = derivatives_path; % Update data paths because I renamed the main folder
[~, ALLEEG_lemon_ec] = pop_loadstudy('filename', [params_lemon.StudyName '.study'], 'filepath', derivatives_path);
etc_lemon_ec = {ALLEEG_lemon_ec.etc};
nRec_lemon_ec = length(ALLEEG_lemon_ec);
clear ALLEEG_lemon_ec

% LEMON eyes open dataset
study = 'LEMON-open';
derivatives_path = '/rechenmagd3/Experiments/2021_preprocessing/datasets/LEMON-8min-open-bids/derivatives_v2023_05_05';
[STUDY_lemon_eo, ALLEEG_lemon_eo] = pop_loadstudy('filename', [study '.study'], 'filepath', derivatives_path);
etc_lemon_eo = {ALLEEG_lemon_eo.etc};
nRec_lemon_eo = length(ALLEEG_lemon_eo);
clear ALLEEG_lemon_eo


% TD-BRAIN dataset
study = 'VanDijk';
derivatives_path = '/rechenmagd3/Experiments/2021_preprocessing/datasets/vanDijk/derivatives_v2023_07_14';
params_vd = define_params(fullfile(derivatives_path,'params.json'));
[STUDY_vd, ALLEEG_vd] = pop_loadstudy('filename', [study '.study'], 'filepath', derivatives_path);

% TD-BRAIN participants.tsv
participants_vd = readtable('/rechenmagd3/Experiments/2021_preprocessing/datasets/vanDijk/participants.tsv','FileType','text','Delimiter', '\t');

% Select only one recording session per participant
filename = {ALLEEG_vd.filename};
ec = find(contains(filename,'restEC'));
filename_ec = filename(contains(filename,'restEC'));
ids_ec = cellfun(@(x) x(1:12),filename_ec,'UniformOutput',false);
[~,ic] = unique(ids_ec);
unique_ec = ec(ic);

eo = find(contains(filename,'restEO'));
filename_eo = filename(contains(filename,'restEO'));
ids_eo = cellfun(@(x) x(1:12),filename_eo,'UniformOutput',false);
[~,io] = unique(ids_eo);
unique_eo = eo(io);

etc_vd_ec = {ALLEEG_vd(unique_ec).etc};
etc_vd_eo = {ALLEEG_vd(unique_eo).etc};
nRec_vd_ec = length(unique_ec);
nRec_vd_eo = length(unique_eo);

clear ALLEEG_vd
%% Study demographics
edges = 0:5:100;
f = figure('Units','centimeters','Position', [0 0 5 8]);
tiledlayout(2,1);

% Age LEMON (N=211)
age_group_lemon = categorical({STUDY_lemon_eo.datasetinfo.age_group});
age_lemon = nan(size(age_group_lemon));
age_lemon(age_group_lemon == '20-25') = 20;
age_lemon(age_group_lemon == '25-30') = 25;
age_lemon(age_group_lemon == '30-35') = 30;
age_lemon(age_group_lemon == '35-40') = 35;
age_lemon(age_group_lemon == '55-60') = 55;
age_lemon(age_group_lemon == '60-65') = 60;
age_lemon(age_group_lemon == '65-70') = 65;
age_lemon(age_group_lemon == '70-75') = 70;
age_lemon(age_group_lemon == '75-80') = 75;

nexttile
histogram(age_lemon,edges,'FaceColor',[0.5 0.5 0.5]);
box off

% Age Van Dijk (N=1256) 18 participants NaN values
age_vd = cellfun(@str2num,participants_vd.age,'UniformOutput',false);
age_vd = cellfun(@(x) x(1),age_vd);

nexttile
histogram(age_vd,edges,'FaceColor',[0.5 0.5 0.5]);
box off
xlabel('Age','FontName','Arial')

% Gender LEMON (N=211)
% 1 = female; 2 = male
gender_lemon = [STUDY_lemon_eo.datasetinfo.sex];
Nfemales_lemon = sum(gender_lemon==1);
Nmales_lemon = sum(gender_lemon==2);

% Gender VanDijk (N=1273) 1 participant NaN value
% 0 = female, 1 = male
gender_vd = participants_vd.gender;
Nfemales_vd = sum(gender_vd==0);
Nmales_vd = sum(gender_vd==1);

saveas(f,fullfile(figures_path,'dataset_demographics.svg'));

%% Recordings layout
f = figure('Units','centimeters','Position', [0 0 5 8]);
tcl = tiledlayout(2,1);
% Lemon dataset layout
nexttile
data = load_preprocessed_data(params_lemon,'sub-010321');
layout = ft_prepare_layout([], data);
ft_plot_layout(layout, 'point', true, 'box', 'no', 'label', false, 'mask', [], 'outline', true, 'pointsymbol','.','pointcolor','k','pointsize',5);
nchans_lemon = size(data.label,1);

% Van Dijk dataset layout
nexttile
data = load_preprocessed_data(params_vd,'sub-19681349_ses-1_task-restEC');
layout = ft_prepare_layout([], data);
ft_plot_layout(layout, 'point', true, 'box', 'no', 'label', false, 'mask', [], 'outline', true, 'pointsymbol','.','pointcolor','k','pointsize',5);
nchans_vd = size(data.label,1);
set(gca,'FontName','Arial');
saveas(f,fullfile(figures_path,'layout_lemon_vanDijk.svg'));
%% Preprocessing summary

% Percentage rejected channels LEMON closed
badchans_lemon_ec = zeros(1,nRec_lemon_ec);
recmask = cellfun(@(x) isfield(x, 'clean_channel_mask'), etc_lemon_ec); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_channel_mask, etc_lemon_ec(recmask),'uni',0);
tmp = cellfun(@(x) x(:), tmp, 'UniformOutput', false); % force to be a column vector
badrecs = double(~cat(2,tmp{:}));
badchans_lemon_ec(recmask) = sum(badrecs,1);
% Average and std
stats.lemon.ec.badchans.mean = mean(badchans_lemon_ec);
stats.lemon.ec.badchans.std = std(badchans_lemon_ec);
badchans_lemon_ec = 100*(badchans_lemon_ec/nchans_lemon);
stats.lemon.ec.badchans.mean_per = mean(badchans_lemon_ec);
stats.lemon.ec.badchans.std_per = std(badchans_lemon_ec);

% Percentage rejected channels LEMON open
badchans_lemon_eo = zeros(1,nRec_lemon_eo);
recmask = cellfun(@(x) isfield(x, 'clean_channel_mask'), etc_lemon_eo); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_channel_mask, etc_lemon_eo(recmask),'uni',0);
tmp = cellfun(@(x) x(:), tmp, 'UniformOutput', false); % force to be a column vector
badrecs = double(~cat(2,tmp{:}));
badchans_lemon_eo(recmask) = sum(badrecs,1);
stats.lemon.eo.badchans.mean = mean(badchans_lemon_eo);
stats.lemon.eo.badchans.std = std(badchans_lemon_eo);
badchans_lemon_eo = 100*(badchans_lemon_eo/nchans_lemon);
stats.lemon.eo.badchans.mean_per = mean(badchans_lemon_eo);
stats.lemon.eo.badchans.std_per = std(badchans_lemon_eo);
badchans_lemon = {badchans_lemon_ec,badchans_lemon_eo};

% Percentage rejected ICs LEMON closed
extract_bad_ICs = @(x) sum((x(:,2:end) > params_lemon.ICLabel(2:end,1)').*(x(:,2:end) < params_lemon.ICLabel(2:end,2)'));
classifications = cellfun(@(x) x.ic_classification.ICLabel.orig_classifications, etc_lemon_ec,'uni',0);
badics = cellfun(extract_bad_ICs, classifications,'uni',0);
classes = etc_lemon_ec{1}.ic_classification.ICLabel.classes;
badics = reshape(cell2mat(badics),length(classes)-1,nRec_lemon_ec);
badics_lemon_ec = sum(badics,1);
stats.lemon.ec.badics.mean = mean(badics_lemon_ec);
stats.lemon.ec.badics.std = std(badics_lemon_ec);
badics_lemon_ec = 100*(badics_lemon_ec./(nchans_lemon-badchans_lemon_ec));
stats.lemon.ec.badics.mean_per = mean(badics_lemon_ec);
stats.lemon.ec.badics.std_per = std(badics_lemon_ec);

% Percentage rejected ICs LEMON open
classifications = cellfun(@(x) x.ic_classification.ICLabel.orig_classifications, etc_lemon_eo,'uni',0);
badics = cellfun(extract_bad_ICs, classifications,'uni',0);
classes = etc_lemon_eo{1}.ic_classification.ICLabel.classes;
badics = reshape(cell2mat(badics),length(classes)-1,nRec_lemon_eo);
badics_lemon_eo = sum(badics,1);
stats.lemon.eo.badics.mean = mean(badics_lemon_eo);
stats.lemon.eo.badics.std = std(badics_lemon_eo);
badics_lemon_eo = 100*(badics_lemon_eo./(nchans_lemon-badchans_lemon_eo));
stats.lemon.eo.badics.mean_per = mean(badics_lemon_eo);
stats.lemon.eo.badics.std_per = std(badics_lemon_eo);
badics_lemon = {badics_lemon_ec,badics_lemon_eo};

% Percentage of rejected time segmets LEMON closed
recmask = cellfun(@(x) isfield(x, 'clean_sample_mask'), etc_lemon_ec); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_sample_mask, etc_lemon_ec(recmask),'uni',0);
badsegs_lemon_ec = cellfun(@(x) sum(~x), tmp);
stats.lemon.ec.badsegs.mean = mean(badsegs_lemon_ec)/250;
stats.lemon.ec.badsegs.std= std(badsegs_lemon_ec)/250;
badsegs_lemon_ec = cellfun(@(x) sum(~x)./length(x), tmp)*100;
stats.lemon.ec.badsegs.mean_per = mean(badsegs_lemon_ec);
stats.lemon.ec.badsegs.std_per = std(badsegs_lemon_ec);

% Percentage of rejected time segmets LEMON open
recmask = cellfun(@(x) isfield(x, 'clean_sample_mask'), etc_lemon_eo); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_sample_mask, etc_lemon_eo(recmask),'uni',0);
badsegs_lemon_eo = cellfun(@(x) sum(~x), tmp);
stats.lemon.eo.badsegs.mean = mean(badsegs_lemon_eo)/250;
stats.lemon.eo.badsegs.std = std(badsegs_lemon_eo)/250;
badsegs_lemon_eo = cellfun(@(x) sum(~x)./length(x), tmp)*100;
stats.lemon.eo.badsegs.mean_per = mean(badsegs_lemon_eo);
stats.lemon.eo.badsegs.std_per = std(badsegs_lemon_eo);
badsegs_lemon = {badsegs_lemon_ec, badsegs_lemon_eo};

% Percentage rejected channels Van Dijk closed
badchans_vd_ec = zeros(1,nRec_vd_ec);
recmask = cellfun(@(x) isfield(x, 'clean_channel_mask'), etc_vd_ec); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_channel_mask, etc_vd_ec(recmask),'uni',0);
tmp = cellfun(@(x) x(:), tmp, 'UniformOutput', false); % force to be a column vector
badrecs = double(~cat(2,tmp{:}));
badchans_vd_ec(recmask) = sum(badrecs,1);
stats.vd.ec.badchans.mean = mean(badchans_vd_ec);
stats.vd.ec.badchans.std = std(badchans_vd_ec);
badchans_vd_ec = 100*(badchans_vd_ec/nchans_vd);
stats.vd.ec.badchans.mean_per = mean(badchans_vd_ec);
stats.vd.ec.badchans.std_per = std(badchans_vd_ec);

% Percentage rejected channels Van Dijk open
badchans_vd_eo = zeros(1,nRec_vd_eo);
recmask = cellfun(@(x) isfield(x, 'clean_channel_mask'), etc_vd_eo); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_channel_mask, etc_vd_eo(recmask),'uni',0);
tmp = cellfun(@(x) x(:), tmp, 'UniformOutput', false); % force to be a column vector
badrecs = double(~cat(2,tmp{:}));
badchans_vd_eo(recmask) = sum(badrecs,1);
stats.vd.eo.badchans.mean = mean(badchans_vd_eo);
stats.vd.eo.badchans.std = std(badchans_vd_eo);
badchans_vd_eo = 100*(badchans_vd_eo/nchans_vd);
stats.vd.eo.badchans.mean_per = mean(badchans_vd_eo);
stats.vd.eo.badchans.std_per = std(badchans_vd_eo);
badchans_vd = {badchans_vd_ec, badchans_vd_eo};

% Percentage rejected ICs VanDijk closed
extract_bad_ICs = @(x) sum((x(:,2:end) > params_vd.ICLabel(2:end,1)').*(x(:,2:end) < params_vd.ICLabel(2:end,2)'));
classifications = cellfun(@(x) x.ic_classification.ICLabel.orig_classifications, etc_vd_ec,'uni',0);
badics = cellfun(extract_bad_ICs, classifications,'uni',0);
classes = etc_vd_ec{1}.ic_classification.ICLabel.classes;
badics = reshape(cell2mat(badics),length(classes)-1,nRec_vd_ec);
badics_vd_ec = sum(badics,1);
stats.vd.ec.badics.mean = mean(badics_vd_ec);
stats.vd.ec.badics.std = std(badics_vd_ec);
badics_vd_ec = 100*(badics_vd_ec./(nchans_vd-badchans_vd_ec));
stats.vd.ec.badics.mean_per = mean(badics_vd_ec);
stats.vd.ec.badics.std_per = std(badics_vd_ec);

% Percentage rejected ICs LEMON open
classifications = cellfun(@(x) x.ic_classification.ICLabel.orig_classifications, etc_vd_eo,'uni',0);
badics = cellfun(extract_bad_ICs, classifications,'uni',0);
classes = etc_vd_eo{1}.ic_classification.ICLabel.classes;
badics = reshape(cell2mat(badics),length(classes)-1,nRec_vd_eo);
badics_vd_eo = sum(badics,1);
stats.vd.eo.badics.mean = mean(badics_vd_eo);
stats.vd.eo.badics.std = std(badics_vd_eo);
badics_vd_eo = 100*(badics_vd_eo./(nchans_vd-badchans_vd_eo));
stats.vd.eo.badics.mean_per = mean(badics_vd_eo);
stats.vd.eo.badics.std_per = std(badics_vd_eo);
badics_vd = {badics_vd_ec,badics_vd_eo};

% Percentage of rejected time segmets Van Dijk closed
recmask = cellfun(@(x) isfield(x, 'clean_sample_mask'), etc_vd_ec); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_sample_mask, etc_vd_ec(recmask),'uni',0);
badsegs_vd_ec = cellfun(@(x) sum(~x), tmp);
stats.vd.ec.badsegs.mean = mean(badsegs_vd_ec)/250;
stats.vd.ec.badsegs.std= std(badsegs_vd_ec)/250;
badsegs_vd_ec = cellfun(@(x) sum(~x)./length(x), tmp)*100;
stats.vd.ec.badsegs.mean_per = mean(badsegs_vd_ec);
stats.vd.ec.badsegs.std_per = std(badsegs_vd_ec);

% Percentage of rejected time segmets Van Dijk open
recmask = cellfun(@(x) isfield(x, 'clean_sample_mask'), etc_vd_eo); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_sample_mask, etc_vd_eo(recmask),'uni',0);
badsegs_vd_eo = cellfun(@(x) sum(~x), tmp);
stats.vd.eo.badsegs.mean = mean(badsegs_vd_eo)/250;
stats.vd.eo.badsegs.std= std(badsegs_vd_eo)/250;
badsegs_vd_eo = cellfun(@(x) sum(~x)./length(x), tmp)*100;
stats.vd.eo.badsegs.mean_per = mean(badsegs_vd_eo);
stats.vd.eo.badsegs.std_per = std(badsegs_vd_eo);
badsegs_vd = {badsegs_vd_ec, badsegs_vd_eo};


% Rainclouds
colors = [0 0.447 0.741; 0.455 0.674 0.188];

f = figure('Units','centimeters','Position', [0 0 12 10]);
tlc = tiledlayout(2,3,'Padding','compact');

nexttile(1);
rm_raincloud_cg(badchans_lemon,'colours',colors,'box_dodge',false,...
    'line_width',1,'raindrop_size',10,'opacity',0.3,'aligned_plots',false,'dist_plots',0.1,'bandwidth',3);
xlim([0 50])
title('Rejected channels (%)');

nexttile(2);
rm_raincloud_cg(badics_lemon,'colours',colors,'box_dodge',false,...
    'line_width',1,'raindrop_size',10,'opacity',0.3,'aligned_plots',false,'dist_plots',0.1,'bandwidth',3);
xlim([0 80])
title('Rejected ICs (%)');

nexttile(3);
rm_raincloud_cg(badsegs_lemon,'colours',colors,'box_dodge',false,...
    'line_width',1,'raindrop_size',10,'opacity',0.3,'aligned_plots',false,'dist_plots',0.1,'bandwidth',3);
xlim([0 80])
title('Rejected time segments (%)');

nexttile(4);
rm_raincloud_cg(badchans_vd,'colours',colors,'box_dodge',false,...
    'line_width',1,'raindrop_size',10,'opacity',0.3,'aligned_plots',false,'dist_plots',0.1,'bandwidth',3);
xlim([0 50])

nexttile(5)
rm_raincloud_cg(badics_vd,'colours',colors,'box_dodge',false,...
    'line_width',1,'raindrop_size',10,'opacity',0.3,'aligned_plots',false,'dist_plots',0.1,'bandwidth',3);
xlim([0 80])

nexttile(6)
rm_raincloud_cg(badsegs_vd,'colours',colors,'box_dodge',false,...
    'line_width',1,'raindrop_size',10,'opacity',0.3,'aligned_plots',false,'dist_plots',0.1,'bandwidth',3);
xlim([0 80])

set(gca,'FontName','Arial');
saveas(f,fullfile(figures_path,'preprocessing_lemon_vanDijk.svg'));
save(fullfile(results_path,'stats_preprocessing_lemon_vanDijk'),'stats');

