%% Datasets demographics and preprocessing.
%
% Code to generate figure 4 of the paper.
%
% Cristina Gil, 10.01.2023, Technical University of Munich


% Load preprocessed data
% params = define_params();
% Add toolboxes to the path
run('/rechenmagd4/toolboxes_and_functions/eeglab/eeglab.m');
addpath('/rechenmagd4/toolboxes_and_functions/fieldtrip');
ft_defaults   
addpath /rechenmagd4/toolboxes_and_functions/plotting_functions/raincloudplots
%% Datasets
% Chronic pain dataset
params.study = 'CP-CGX';
params.raw_data_path = '/rechenmagd4/Experiments/2022_chronic_pain_CGX/cpCGX_BIDS';
t = '2022_12_12';
params.preprocessed_data_path = fullfile(params.raw_data_path, ['derivatives_v' t]);
[STUDY, ALLEEG] = pop_loadstudy('filename', [params.study '.study'], 'filepath', params.preprocessed_data_path);
etc_cp = {ALLEEG.etc};
nRec_cp = length(ALLEEG);
age_cp = str2double({STUDY.datasetinfo.age});
sex_cp =  categorical({STUDY.datasetinfo.sex});
clear ALLEEG STUDY;

% LEMON dataset
params.study = 'LEMON-2s';
params.raw_data_path = '/rechenmagd4/Experiments/2021_preprocessing/datasets/LEMON-8min-bids';
t = '2022_12_14';
params.preprocessed_data_path = fullfile(params.raw_data_path, ['derivatives_v' t]);
[STUDY_lemon, ALLEEG_lemon] = pop_loadstudy('filename', [params.study '.study'], 'filepath', params.preprocessed_data_path);
etc_lemon = {ALLEEG_lemon.etc};
nRec_lemon = length(ALLEEG_lemon);
age_group = categorical({STUDY_lemon.datasetinfo.age_group});
subject = categorical({STUDY_lemon.datasetinfo.subject});
clear ALLEEG_lemon STUDY_lemon

% TD-BRAIN dataset
study = 'Van-Dijk';
preprocessed_data_path = '/rechenmagd3/Experiments/2021_preprocessing/datasets/vanDijk/derivatives_v2023_07_14';
[STUDY_vd, ALLEEG_vd] = pop_loadstudy('filename', [study '.study'], 'filepath', preprocessed_data_path);
etc_vd = {ALLEEG_vd.etc};
nRec_vd = length(ALLEEG_vd);
age_vd = categorical({STUDY_lemon.datasetinfo.age_group});
subject_vd = categorical({STUDY_lemon.datasetinfo.subject});
%% Study demographics
edges = 0:5:100;
f = figure;
tcl = tiledlayout(2,1);

% Age CP (N = 74)
nexttile
histogram(age_cp,edges)
box off

% Age LEMON (N=210)
age_lemon = nan(size(age_group));
age_lemon(age_group == '20-25') = 20;
age_lemon(age_group == '25-30') = 25;
age_lemon(age_group == '30-35') = 30;
age_lemon(age_group == '35-40') = 35;
age_lemon(age_group == '55-60') = 55;
age_lemon(age_group == '60-65') = 60;
age_lemon(age_group == '65-70') = 65;
age_lemon(age_group == '70-75') = 70;
age_lemon(age_group == '75-80') = 75;

nexttile
histogram(age_lemon,edges);
box off

% Sex CP
summary(sex_cp);

% Sex LEMON (N=212)
demo = readtable('/rechenmagd4/Experiments/2021_preprocessing/datasets/LEMON/Participants_MPILMBB_LEMON.csv');
demo.Properties.VariableNames = {'participant_id','sex','age_group'};
demo.participant_id = categorical(demo.participant_id);
% Add subjects whose data was truncated
subject(end+1) = 'sub-010015';
subject(end+1) = 'sub-010100';
s = demo(ismember(demo.participant_id,subject),:);
age_group = categorical(s.age_group);
sex_cp = categorical(s.sex);
summary(sex_cp);
mask_young = ismember(age_group,{'20-25','25-30','30-35'});
summary(sex_cp(mask_young))
mask_old = ismember(age_group,{'55-60','60-65','65-70','70-75','75-80'});
summary(sex_cp(mask_old));

%% Recordings layout
% Lemon dataset layout
data = load_preprocessed_data(params,'sub-010321');
f = figure;
layout = ft_prepare_layout([], data);
ft_plot_layout(layout, 'point', true, 'box', 'no', 'label', false, 'mask', [], 'outline', true);
saveas(f,'/media/cgil/Verbatim/figures/layout_lemon.svg');

% Chronic pain layout
params.study = 'CP-CGX';
params.raw_data_path = '/rechenmagd4/Experiments/2022_chronic_pain_CGX/cpCGX_BIDS';
t = '2022_12_12';
params.preprocessed_data_path = fullfile(params.raw_data_path, ['derivatives_v' t]);
data = load_preprocessed_data(params,'sub-001_task-EC');
f = figure;
layout = ft_prepare_layout([], data);
ft_plot_layout(layout, 'point', true, 'box', 'no', 'label', false, 'mask', [], 'outline', true);
saveas(f,'/media/cgil/Verbatim/figures/layout_CGX.svg');
%% Preprocessing
% Nr rejected channels CP
badchans_cp = zeros(1,nRec_cp);
recmask = cellfun(@(x) isfield(x, 'clean_channel_mask'), etc_cp); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_channel_mask, etc_cp(recmask),'uni',0);
tmp = cellfun(@(x) x(:), tmp, 'UniformOutput', false); % force to be a column vector
badrecs = double(~cat(2,tmp{:}));
badchans_cp(recmask) = sum(badrecs,1);
badchans_cp = 100*(badchans_cp/29);

% Nr rejected channels LEMON
badchans_lemon = zeros(1,nRec_lemon);
recmask = cellfun(@(x) isfield(x, 'clean_channel_mask'), etc_lemon); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_channel_mask, etc_lemon(recmask),'uni',0);
tmp = cellfun(@(x) x(:), tmp, 'UniformOutput', false); % force to be a column vector
badrecs = double(~cat(2,tmp{:}));
badchans_lemon(recmask) = sum(badrecs,1);
badchans_lemon = 100*(badchans_lemon/61);

% Nr rejected ICs CP
extract_bad_ICs = @(x) sum((x(:,2:end) > params.IClabel(2:end,1)').*(x(:,2:end) < params.IClabel(2:end,2)'));
classifications = cellfun(@(x) x.ic_classification.ICLabel.orig_classifications, etc_cp,'uni',0);
badics = cellfun(extract_bad_ICs, classifications,'uni',0);
classes = etc_cp{1}.ic_classification.ICLabel.classes;
badics = reshape(cell2mat(badics),length(classes)-1,nRec_cp);
badics_cp = sum(badics,1);
badics_cp = 100*(badics_cp./(29-badchans_cp)); % Percentage of rejected ICs out of all clean channels

% Nr rejected ICs LEMON
classifications = cellfun(@(x) x.ic_classification.ICLabel.orig_classifications, etc_lemon,'uni',0);
badics = cellfun(extract_bad_ICs, classifications,'uni',0);
classes = etc_lemon{1}.ic_classification.ICLabel.classes;
badics = reshape(cell2mat(badics),length(classes)-1,nRec_lemon);
badics_lemon = sum(badics,1);
badics_lemon = 100*(badics_lemon./(61-badchans_lemon));

% Percentage of rejected time segments CP
recmask = cellfun(@(x) isfield(x, 'clean_sample_mask'), etc_cp); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_sample_mask, etc_cp(recmask),'uni',0);
badsegs_cp = cellfun(@(x) sum(~x)./length(x), tmp)*100;

% Percentage of rejected time segmets LEMON
recmask = cellfun(@(x) isfield(x, 'clean_sample_mask'), etc_lemon); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_sample_mask, etc_lemon(recmask),'uni',0);
badsegs_lemon = cellfun(@(x) sum(~x)./length(x), tmp)*100;

% Rainclouds
colors = lines(2);
f = figure('Units','centimeters','Position', [0 0 10 20]);
tlc = tiledlayout(2,3,'Padding','compact');

nexttile(1);
rm_raincloud_cg({badchans_cp},'colours',colors(1,:),'bxcl',colors(1,:),'box_dodge',false,...
    'line_width',1,'raindrop_size',20,'opacity',0.4,'aligned_plots',false,'dist_plots',0.1,'bandwidth',3);
xlim([0 50])
title('Percentage of rejected channels');

nexttile(4);
rm_raincloud_cg({badchans_lemon},'colours',colors(1,:),'bxcl',colors(1,:),'box_dodge',false,...
    'line_width',1,'raindrop_size',20,'opacity',0.4,'aligned_plots',false,'dist_plots',0.1,'bandwidth',3);
xlim([0 50])

nexttile(2);
rm_raincloud_cg({badics_cp},'colours',colors(1,:),'bxcl',colors(1,:),'box_dodge',false,...
    'line_width',1,'raindrop_size',20,'opacity',0.4,'aligned_plots',false,'dist_plots',0.1,'bandwidth',3);
xlim([0 50])
title('Percentage of rejected ICs');

nexttile(5)
rm_raincloud_cg({badics_lemon},'colours',colors(1,:),'bxcl',colors(1,:),'box_dodge',false,...
    'line_width',1,'raindrop_size',20,'opacity',0.4,'aligned_plots',false,'dist_plots',0.1,'bandwidth',3);
xlim([0 50])

nexttile(3);
rm_raincloud_cg({badsegs_cp},'colours',colors(1,:),'bxcl',colors(1,:),'box_dodge',false,...
    'line_width',1,'raindrop_size',20,'opacity',0.4,'aligned_plots',false,'dist_plots',0.1,'bandwidth',3);
xlim([0 80])
title('Percentage of rejected time segments');

nexttile(6)
rm_raincloud_cg({badsegs_lemon},'colours',colors(1,:),'bxcl',colors(1,:),'box_dodge',false,...
    'line_width',1,'raindrop_size',20,'opacity',0.4,'aligned_plots',false,'dist_plots',0.1,'bandwidth',3);
xlim([0 80])

saveas(f,'/media/cgil/Verbatim/figures/preprocessing_rc.svg');