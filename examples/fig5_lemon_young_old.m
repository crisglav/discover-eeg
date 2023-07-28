%% Comparison of EEG features between old and young groups on the LEMON dataset
%
% This script compares the EEG features obtained after applyting the
% pipeline between old and young people on the LEMON dataset eyes closed 
% and generates Figure 5 of the manuscript.
%
% The bayesFactor package for Bayesian statistics testing 
% and the Raincloud plots function need to be downloaded for statistical
% analysis and visualization, respectively. The code can be found in:
%
% - BayesFactor [https:/www.github.com/klabhub/bayesFactor]
% - Raincloud plots [https:/www.github.com/crisglav/Raincloud_RM]
%
% Cristina Gil, 21.11.2022, Technical University of Munich

clear all;

% Add toolboxes and functions
eeglab_path = '/rechenmagd4/toolboxes_and_functions/eeglab';
run(fullfile(eeglab_path,'eeglab.m'));
addpath /rechenmagd4/toolboxes_and_functions/fieldtrip
ft_defaults
addpath(genpath('../external_functions'));
run('installBayesFactor.m')

figures_path = '/rechenmagd3/Experiments/2021_preprocessing/figures';
results_path = '/rechenmagd3/Experiments/2021_preprocessing/results';

% Load atlas for plotting the connectivity matrices
n_sources = 100;
atlas_path = '../parcellations/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv';
atlas100 = readtable(atlas_path);

%% Load study with clean data
study_path = '/rechenmagd3/Experiments/2021_preprocessing/datasets/LEMON-8min-closed-bids/derivatives_v2022_12_14';
STUDY = pop_loadstudy('filename', 'LEMON-2s-clean.study', 'filepath', study_path);

% Sort participants by age group
subject = {STUDY.datasetinfo.subject};
age_group = categorical({STUDY.datasetinfo.age_group});
summary(age_group)

% young (20-35)
mask_young = ismember(age_group,{'20-25','25-30','30-35'});
n_young = sum(mask_young);
young_ids = subject(mask_young);

% old (55-80)
mask_old = ismember(age_group,{'55-60','60-65','65-70','70-75','75-80'});
n_old = sum(mask_old);
old_ids = subject(mask_old);

freqBands = {'theta','alpha','beta','gamma'};
connMeas = {'dwpli','aec'};
networks = {'Vis','SomMot','DorsAttn','SalVentAttn','Limbic','Cont','Default'};
% For plotting connectivity matrices
pos = cell(1,length(networks));
axisticks = ones(1,length(networks)+1);
axistickslabelpos = ones(1,length(networks));
for i=1:length(networks)
    pos{i}  = find(cellfun(@(x) contains(x,['_' networks{i} '_']), atlas100.ROIName));
    % This is for plotting the separation lines between networks in the
    % connectivity matrices
    axisticks(i+1) = length(pos{i})+axisticks(i);
    axistickslabelpos(i) = axisticks(i)+(axisticks(i+1)-axisticks(i))/2;
end
newpos = vertcat(pos{:});

%% ALPHA PEAK FREQUENCY
% Load peak frequency data from all recordings
apf_path = fullfile(study_path, 'EEG_features','power');

nSubj = length(subject);
apf_max = nan(1,nSubj);
apf_cog = nan(1,nSubj);

for i=1:nSubj
    load(fullfile(apf_path,[subject{i} '_peakfrequency.mat']));
    if ~isempty(peakfrequency.localmax)
        apf_max(i) = peakfrequency.localmax;
    end
    apf_cog(i) = peakfrequency.cog;
end

% Bayesian one-side independent samples t-test to check if the old group old has lower APF than the young group
young_max = apf_max(mask_young);
old_max = apf_max(mask_old);
[bf_max,p_localmax] = bf.ttest2(old_max,young_max);

young_cog = apf_cog(mask_young);
old_cog = apf_cog(mask_old);
[bf_cog,p_cog] = bf.ttest2(old_cog,young_cog);
save(fullfile(results_path,'stats_lemon_young_old_apf.mat'),'bf_max','p_localmax','bf_cog','p_cog');

%% APF - raincloudplots
colours = lines(2);
f = figure('Units','centimeters','Position', [0 0 16 7]);
tlc = tiledlayout(1,2,'Padding','compact');
ax = nexttile;
rc_max = rm_raincloud({young_max,old_max},colours);
xlim([7.5, 13.5]);
ylim([-0.4,1]);
xlabel('Frequency (Hz)')
title('Local maximum');
nexttile
rc_cog = rm_raincloud({young_cog,old_cog},colours);
xlim([7.5, 13.5]);
ylim([-0.4,1]);
xlabel('Frequency (Hz)')
title('Center of gravity');
hold on
% fake legend
h(1) = scatter(nan,nan,10,colours(1,:),'filled');
hold on
h(2) = scatter(nan,nan,10,colours(2,:),'filled');
legend(h,'young','old');
saveas(f,fullfile(figures_path,'lemon_young_old_apf.svg'));

%% CONNECTIVITY
% Path to connectivity matrices
connectivity_path = fullfile(study_path, 'EEG_features','connectivity');

for iConnMeas = 1: length(connMeas)
    meas = connMeas{iConnMeas};
    % List matrices
    conn_files = dir(fullfile(connectivity_path,['*_' meas '_*.mat']));
    for iBand=1:length(freqBands)
        fBand = freqBands{iBand};
        conn_band = conn_files(contains({conn_files.name},fBand));
        subject_conn = cellfun(@(x) regexp(x,['.*(?=_' meas '_*)'],'match','lineanchors'),{conn_band.name});
        
        conn_old_ids = find(contains(subject_conn,old_ids));
        conn_young_ids = find(contains(subject_conn,young_ids));
        
        old_conn = nan(n_sources,n_sources,length(conn_old_ids));
        young_conn = nan(n_sources,n_sources,length(conn_young_ids));
        
        % Load connectivity matrices from the old group
        for i=1:length(conn_old_ids)
            try
                load(fullfile(conn_band(conn_old_ids(i)).folder,conn_band(conn_old_ids(i)).name));
            catch
                error(['cannot load data from ' conn_band(conn_old_ids(i)).name])
            end
            old_conn(:,:,i) = connMatrix;
        end
        % Load connectivity matrices from the young group
        for i=1:length(conn_young_ids)
            try
                load(fullfile(conn_band(conn_young_ids(i)).folder,conn_band(conn_young_ids(i)).name));
            catch
                error(['cannot load data from ' conn_band(conn_young_ids(i)).name])
            end
            young_conn(:,:,i) = connMatrix;
        end
        
        % Re-arrange matrices by atlas networks
        old_conn_r = old_conn(newpos,newpos,:);
        young_conn_r = young_conn(newpos,newpos,:);
        
        % 9900 bayesian two-sided independent sample ttests between old and young
        % connectivity matrices (one test per source pair)
        b = nan(n_sources,n_sources);
        pval = nan(n_sources,n_sources);
        tval = nan(n_sources,n_sources);
        for i=2:n_sources
            for j=1:(i-1)
                o = squeeze(old_conn_r(i,j,:));
                y = squeeze(young_conn_r(i,j,:));
                [~,p,~,s] = ttest2(o,y);
                pval(i,j) = p;
                tval(i,j) = s.tstat;
                b(i,j) = bf.ttest('T',s.tstat,'N',[numel(o) numel(y)]);
                
            end
        end
        
        % Corret p-values for multiple comparisons with fdr
        [~,~,padj] = fdr(pval(:));
        padj = reshape(padj,n_sources,n_sources);
        
        stats.(fBand).tval = tval;
        stats.(fBand).pval = pval;
        stats.(fBand).padj = padj;
        stats.(fBand).bf = b;
        
    end
    save(fullfile(results_path,['stats_lemon_young_old_conn_' meas '.mat']),'stats');
end
%% Connectivity - Plot t-values and Bayes factors
cmap = turbo;
cmap = cmap(20:235,:); % limit the colormap so that it does not reach very dark red and very dark blue
for iConnMeas = 1: length(connMeas)
    % Load stats file
    meas = connMeas{iConnMeas};
    statsfile = fullfile(results_path,['stats_lemon_young_old_conn_' meas '.mat']);
    load(statsfile);
    
    fig = figure('Units','centimeters','Position',[0 0 10 10], 'visible', 'on');
    tcl = tiledlayout(2,2,'TileSpacing','compact','Padding','none');
    cmin = -4;
    cmax = 4;
    % Plot matrices with t-values color coded. Shade in grey all t-values
    % with low evidence (1/30 < BF < 30)
    for iBand=1:length(freqBands)
        fBand = freqBands{iBand};
        
        % Bayesian statistics
        high_evidence_cc = or(stats.(fBand).bf >= 30, stats.(fBand).bf <= 1/30);
        low_evidence = and(stats.(fBand).bf > 1/30, stats.(fBand).bf < 30);
        
        %     % Frequentist statistics
        %     high_evidence = stats.(fBand).padj < 0.05;
        %     low_evidence = stats.(fBand).padj >= 0.05;
        
        ax = nexttile(tcl);
        imagesc(stats.(fBand).tval,'AlphaData',0.2,[cmin cmax]);
        colormap(ax,cmap);
        hold on
        imagesc(stats.(fBand).tval,'AlphaData',high_evidence_cc, [cmin cmax]);
        hold off
        % Set to white the upper triangular matrix
        ax2 = copyobj(ax,ax.Parent);
        im = imagesc(ax2,ones(size(stats.(fBand).tval))*100,'AlphaData',triu(ones(size(stats.(fBand).tval)))); 
        colormap(ax2,white);
        set(ax2,'Color','none');
        % Plot grid and network labels
        set(ax,'XTick',axistickslabelpos,'XtickLabel',[],'YTick',axistickslabelpos,'YtickLabel',[],'TickLength',[0 0],'TickDir','out','box','off')
        if iBand ==3
            set(ax,'XtickLabel',networks,'XtickLabelRotation',45,'YtickLabel',networks)
        end
        grid(ax2,'on');
        set(ax2,'GridColor','w','GridAlpha',1,'Layer','top');
        set(ax2,'Ytick',axisticks','Xtick',axisticks,'yticklabel',[],'xticklabel',[]);
        set(ax2,'TickLength',[0 0],'TickDir','out','box','off')
        
%         colorbar;
        title(fBand);
    end
    title(tcl,meas);

    saveas(fig,fullfile(figures_path,['lemon_young_old_conn_' meas '.svg']));
end
% % Generate colorbar
% f = figure('Units','centimeters','Position',[0 0 3 4]);
% colormap(cmap)
% ax = axes;
% c = colorbar(ax,'Location','eastoutside');
% caxis([cmin cmax]);
% ax.Visible = 'off';
% saveas(f,'colorbar.svg');
%% GRAPH MEASURES
graph_path = fullfile(study_path, 'EEG_features','graph_measures');

for iConnMeas = 1: length(connMeas)
    % List matrices
    meas = connMeas{iConnMeas};
    graph_files = dir(fullfile(graph_path,['*_' meas '_*.mat']));
    
    for iBand=1:length(freqBands)
        fBand = freqBands{iBand};
        graph_files_band = graph_files(contains({graph_files.name},fBand));
        subject_graph = cellfun(@(x) regexp(x,['.*(?=_graph_' meas '_*)'],'match','lineanchors'),{graph_files_band.name});
        
        % Load graph measures from the old and young group
        graph_old_ids = find(contains(subject_graph,old_ids));
        graph_young_ids = find(contains(subject_graph,young_ids));
        
        old_degree = nan(length(graph_old_ids),n_sources);
        old_cc = nan(length(graph_old_ids),n_sources);
        old_gcc = nan(1,length(graph_old_ids));
        old_geff = nan(1,length(graph_old_ids));
        old_s = nan(1,length(graph_old_ids));
        young_degree = nan(length(graph_young_ids),n_sources);
        young_cc = nan(length(graph_young_ids),n_sources);
        young_gcc = nan(1,length(graph_young_ids));
        young_geff = nan(1,length(graph_young_ids));
        young_s = nan(1,length(graph_young_ids));
        
        for i=1:length(graph_old_ids)
            try
                load(fullfile(graph_files_band(graph_old_ids(i)).folder,graph_files_band(graph_old_ids(i)).name));
            catch
                error(['cannot load data from ' graph_files_band(graph_old_ids(i)).name])
            end
            old_degree(i,:) = graph_measures.degree;
            old_cc(i,:) = graph_measures.cc;
            old_gcc(i) = graph_measures.gcc;
            old_geff(i) = graph_measures.geff;
            old_s(i) = graph_measures.smallworldness;
        end
        
        for i=1:length(graph_young_ids)
            try
                load(fullfile(graph_files_band(graph_young_ids(i)).folder,graph_files_band(graph_young_ids(i)).name));
            catch
                error(['cannot load data from ' graph_files_band(graph_young_ids(i)).name])
            end
            young_degree(i,:) = graph_measures.degree;
            young_cc(i,:) = graph_measures.cc;
            young_gcc(i) = graph_measures.gcc;
            young_geff(i) = graph_measures.geff;
            young_s(i) = graph_measures.smallworldness;
        end
        
        data.gcc(iBand,:) = {young_gcc,old_gcc};
        data.geff(iBand,:) = {young_geff,old_geff};
        data.s(iBand,:) = {young_s,old_s};
        
        % Statistics - local measures
        t_degree = nan(1,n_sources);
        b_degree = nan(1,n_sources);
        p_degree = nan(1,n_sources);
        
        t_cc = nan(1,n_sources);
        b_cc = nan(1,n_sources);
        p_cc = nan(1,n_sources);
        % 100 two sided independent sample Bayesian ttests between old and young
        % local measures
        for i=1:n_sources
            o_degree = old_degree(:,i);
            y_degree = young_degree(:,i);
            [~,p,~,s] = ttest2(o_degree,y_degree);
            p_degree(i) = p;
            t_degree(i) = s.tstat;
            b_degree(i) = bf.ttest('T',s.tstat,'N',[numel(o_degree) numel(y_degree)]);
            
            o_cc = old_cc(:,i);
            y_cc = young_cc(:,i);
            [~,p,~,s] = ttest2(o_cc,y_cc);
            p_cc(i) = p;
            t_cc(i) = s.tstat;
            b_cc(i) = bf.ttest('T',s.tstat,'N',[numel(o_cc) numel(y_cc)]);
        end
        % Corret p-values for multiple comparisons with fdr
        [~,~,padj_degree] = fdr(p_degree);
        [~,~,padj_cc]     = fdr(p_cc);
        
        stats.(fBand).degree.tval = t_degree;
        stats.(fBand).degree.pval = p_degree;
        stats.(fBand).degree.padj = padj_degree;
        stats.(fBand).degree.bf   = b_degree;
        stats.(fBand).cc.tval = t_cc;
        stats.(fBand).cc.pval = p_cc;
        stats.(fBand).cc.padj = padj_cc;
        stats.(fBand).cc.bf   = b_cc;
        
        % Statistics - global measures
        % Two-sided independent samples Bayesian t-tests between old and young global measures
        
        [~,p,~,s] = ttest2(old_gcc,young_gcc);
        b = bf.ttest('T',s.tstat,'N',[numel(old_gcc) numel(young_gcc)]);
        stats.(fBand).gcc.tval = s.tstat;
        stats.(fBand).gcc.pval = p;
        stats.(fBand).gcc.bf = b;
        
        [~,p,~,s] = ttest2(old_geff,young_geff);
        b = bf.ttest('T',s.tstat,'N',[numel(old_geff) numel(young_geff)]);
        stats.(fBand).geff.tval = s.tstat;
        stats.(fBand).geff.pval = p;
        stats.(fBand).geff.bf = b;
        
        [~,p,~,s] = ttest2(old_s,young_s);
        b = bf.ttest('T',s.tstat,'N',[numel(old_s) numel(young_s)]);
        stats.(fBand).s.tval = s.tstat;
        stats.(fBand).s.pval = p;
        stats.(fBand).s.bf = b;
        
    end
    save(fullfile(results_path,['stats_lemon_young_old_graph_' meas '.mat']),'stats');
    save(fullfile(results_path,['data_lemon_young_old_graph_' meas '.mat']),'data');
end
%% Raincloud plots - global measures
colours = repmat(lines(2), [1 1 4]);

for iConnMeas = 1: length(connMeas)
    meas = connMeas{iConnMeas};
    datafile = fullfile(results_path,['data_lemon_young_old_graph_' meas '.mat']);
    load(datafile);
    
    f = figure('Units','centimeters','Position', [0 0 18 5]);
    tlc = tiledlayout(1,3,'Padding','compact');
    ax = nexttile;
    rm_raincloud_cg(data.gcc,'colours',colours,'bandwidth',[],'plot_median_lines',false,...
        'line_width',1,'raindrop_size',10,'opacity',0.4);
    box off
    ax.Color = 'none';
    title('Global cc');
    set(ax,'YTickLabel',flip(freqBands));
    ax = nexttile;
    rm_raincloud_cg(data.geff,'colours',colours,'bandwidth',[],'plot_median_lines',false,...
        'line_width',1,'raindrop_size',10,'opacity',0.4);
    title('Global efficiency');
    set(ax,'YTickLabel',flip(freqBands));
    box off
    ax.Color = 'none';
    ax = nexttile;
    rm_raincloud_cg(data.s,'colours',colours,'bandwidth',[],'plot_median_lines',false,...
        'line_width',1,'raindrop_size',10,'opacity',0.4);
    title('Global smallworldness');
    set(ax,'YTickLabel',flip(freqBands));
    box off
    ax.Color = 'none';
    % fake legend
    h(1) = scatter(nan,nan,10,colours(1,:,1),'filled');
    hold on
    h(2) = scatter(nan,nan,10,colours(2,:,1),'filled');
    legend(h,'young','old');
    saveas(f,fullfile(figures_path,['lemon_young_old_graph_global_ ' meas '.svg']));
end
%% Brain plots - local measures
% Load surface structure
surf = ft_read_headshape('surface_white_both.mat');
% Source model: centroid positons from Schaefer atlas
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = [atlas100.R, atlas100.A, atlas100.S];
cfg.unit = 'mm';
sourcemodel_atlas = ft_prepare_sourcemodel(cfg);
sourcemodel_atlas.coordsys = 'mni';

for iConnMeas = 1: length(connMeas)
    
    meas = connMeas{iConnMeas};
    statsfile = fullfile(results_path,['stats_lemon_young_old_graph_' meas '.mat']);
    load(statsfile);
    
    f_degree = figure('Units','centimeters','Position',[0 0 10 9]);
    f_cc = figure('Units','centimeters','Position',[0 0 10 9]);
    tlc_degree = tiledlayout(f_degree,2,2,'TileSpacing','none','Padding','compact');
    tlc_cc = tiledlayout(f_cc,2,2,'TileSpacing','none','Padding','compact');
    
    for iBand=1:length(freqBands)
        
        fBand = freqBands{iBand};
        cmin = -4;
        cmax = 4;
        
        % Bayesian statistics
        high_evidence_degree = or(stats.(fBand).degree.bf >= 30, stats.(fBand).degree.bf <= 1/30);
        high_evidence_cc = or(stats.(fBand).cc.bf >= 30, stats.(fBand).cc.bf <= 1/30);
        
        %     % Frequentist statistics
        %     high_evidence_degree = stats.(fBand).degree.padj < 0.05;
        %     high_evidence_cc = stats.(fBand).cc.padj < 0.05;
        
        % Plot t-values color coded, only high evidence/significant values
        % Degree
        index = fix((stats.(fBand).degree.tval-cmin)/(cmax-cmin)*length(cmap))+1;
        rgb = squeeze(ind2rgb(index,cmap));
        ax1 = nexttile(tlc_degree);
        ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.2);
        ft_plot_mesh(sourcemodel_atlas.pos(high_evidence_degree,:), 'vertexsize',12, 'vertexcolor',rgb(high_evidence_degree,:));
        ax1.CLim = [cmin cmax];
        colormap(ax1,cmap);
        title(fBand);
       
        % Global clustering coefficient
        index = fix((stats.(fBand).cc.tval-cmin)/(cmax-cmin)*length(cmap))+1;
        rgb = squeeze(ind2rgb(index,cmap));
        ax2 = nexttile(tlc_cc);
        ft_plot_mesh(surf, 'edgecolor', 'none', 'vertexcolor', 'curv','facealpha',0.2);
        ft_plot_mesh(sourcemodel_atlas.pos(high_evidence_cc,:), 'vertexsize',12, 'vertexcolor',rgb(high_evidence_cc,:));
        ax2.CLim = [cmin cmax];
        colormap(ax2,cmap);
        title(fBand);
    end
    colorbar(ax1);
    colorbar(ax2);
    saveas(f_degree,fullfile(figures_path,['lemon_young_old_graph_degree_ ' meas '.bmp']));
    saveas(f_cc,fullfile(figures_path,['lemon_young_old_graph_cc_' meas '.bmp']));
end
%% colorbar figure
% f = figure('Units','centimeters','Position',[0 0 3 4]);
% ax = axes;
% c = colorbar(ax,'Location','southoutside');
% caxis([-4 4]);
% ax.Visible = 'off';
% saveas(f,'colorbar.svg');