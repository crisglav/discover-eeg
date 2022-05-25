%
% 
% Requirements:
% - EEGLab
% - Dipfit (just electrode template)
% - ICLabel
% - CleanRawData
function main_pipeline()

clear all; close all;
rng('default'); % For reproducibility - See discussion in https://sccn.ucsd.edu/pipermail/eeglablist/2022/016932.html

% Define the parameters
params = define_params();
save(fullfile(params.preprocessed_data_path,'pipeline_params.mat'),'params');

%% Import data
% Import raw data in BIDS format...
[STUDY, ALLEEG, bids] = pop_importbids(params.raw_data_path,'outputdir',params.preprocessed_data_path,...
    'studyName',params.study,'bidstask',params.task,'bidschanloc','on','bidsevent','off');

% Select only EEG channels
eegchans = find(contains(lower({ALLEEG(1).chanlocs.type}),'eeg'));
ALLEEG = pop_select(ALLEEG , 'channel', eegchans);

CURRENTSTUDY = 1;
EEG = ALLEEG;
CURRENTSET = 1:length(EEG);

% If there are no channel locations loaded use a default template
if isempty(find(~cellfun('isempty', {ALLEEG(1).chanlocs.X }), 1))
    EEG=pop_chanedit(EEG, 'lookup',params.elec_template);
end

%% ======== PREPROCESSING =========

%% 1. REMOVE BAD CHANNELS
[EEG.urchanlocs] = deal(EEG.chanlocs); % Keep original channels
EEG = pop_clean_rawdata(EEG,'FlatlineCriterion', params.FlatlineCriterion,...
                            'ChannelCriterion',params.ChannelCriterion,...
                            'LineNoiseCriterion',params.LineNoiseCriterion,...
                            'Highpass',params.Highpass,...
                            'FuseChanRej','on',... % I don't know if this functionality works
                            'BurstCriterion','off',...
                            'WindowCriterion','off',...
                            'BurstRejection','off',...
                            'Distance','Euclidian',...
                            'WindowCriterionTolerances','off');

% Visualization of bad channels
plot_badchannels(params,EEG);

%% 2. REREFERENCE TO AVERAGE REFERENCE
EEG = pop_reref(EEG,[],'interpchan',[]);

%% 3. REMOVE ARTIFACTS WITH ICA
EEG = pop_runica(EEG,'icatype','runica','concatcond','on');
EEG = pop_iclabel(EEG,'default');
EEG = pop_icflag(EEG, params.IClabel); % flag artifactual components using IClabel
classifications = cellfun(@(x) x.ic_classification.ICLabel.classifications, {EEG.etc},'UniformOutput', false); % Keep classifications before component substraction
EEG = pop_subcomp(EEG,[],0); % Substract artifactual independent components
for iRec=1:length(classifications)
    EEG(iRec).etc.ic_classification.ICLabel.orig_classifications = classifications{iRec};
end

% Visualization of rejected ICs
plot_ICs(params,EEG);

%% 4. INTERPOLATE MISSING CHANNELS
for iRec=1:size(EEG,2)
    if ~isempty(EEG(iRec).chaninfo.removedchans)
        l = cellfun(@(x) strcmp(x,'average'), {EEG(iRec).chaninfo.removedchans.ref});  % Deal with a small bug that duplicates bad channels after rereferencing in removedchans
        EEG(iRec).chaninfo.removedchans = EEG(iRec).chaninfo.removedchans(~l);
        EEGtemp = eeg_checkset(EEG(iRec),'loaddata'); % Retrieve data
        EEGtemp = pop_interp(EEGtemp, EEGtemp.chaninfo.removedchans, 'spherical');
        EEGtemp = pop_saveset(EEGtemp, 'savemode', 'resave');
        EEGtemp.data = 'in set file'; % clear data from memory
        EEG = eeg_store(EEG, EEGtemp, iRec);
    end
end

%% 5. REMOVE BAD TIME SEGMENTS
EEG = pop_clean_rawdata(EEG,'FlatlineCriterion','off',...
                            'ChannelCriterion','off',...
                            'LineNoiseCriterion','off',...
                            'Highpass','off',...
                            'BurstCriterion',params.BurstCriterion,...
                            'WindowCriterion',params.WindowCriterion,...
                            'BurstRejection','on',...
                            'Distance','Euclidian',...
                            'WindowCriterionTolerances',params.WindowCriterionTolerances);

% Visualization of rejected time segments
plot_badtimesegments(params,EEG);

%% 6. SEGMENT DATA INTO EPOCHS
% Note: Epochs containing discontinutities will be automatically rejected.
for iRec=1:size(EEG,2)
    % Create markers each x seconds and add them at the end of existing markers
    epoch_start = 1 : EEG(iRec).srate * params.epoch_length * (1-params.epoch_overlap) : EEG(iRec).pnts;
    nEpochs = length(epoch_start);
    if ~isempty(EEG(iRec).event)
        events = EEG(iRec).event;
        n = length(events);        
        [events(n+1:n+nEpochs).type] = deal('epoch_start');  
        latency = num2cell([[events.latency],epoch_start]);
        [events.latency] = latency{:};
        [events(n+1:n+nEpochs).duration] = deal(1);
        
    else % Deal with the case in which no bad time segments were detected
        events = struct('type',repmat({'epoch_start'},1,nEpochs),...
                        'latency',num2cell(epoch_start),...
                        'duration',num2cell(ones(1,nEpochs)));
    end

    EEG(iRec).event = events;
end
EEG = pop_epoch(EEG,{'epoch_start'},[0 params.epoch_length],'epochinfo','yes');

% Save study
EEG = eeg_checkset(EEG);
EEG = pop_saveset(EEG,'savemode','resave');
STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', [params.study '_preprocessed'], 'filepath', params.preprocessed_data_path);

clear EEG ALLEEG;
% ======== END OF PREPROCESSING =========

%% ======= EXTRACTION OF BRAIN FEATURES =========

% You can start directly with preprocessed data in BIDS format by loading an EEGLAB STUDY
% [STUDY, EEG] = pop_loadstudy('filename', [params.study '_preprocessed.study'], 'filepath', params.preprocessed_data_path);
%%
for iRec=1:length(STUDY.datasetinfo)
    
    % BIDS ID
    x = strsplit(STUDY.datasetinfo(iRec).filename,'_eeg.set');
    bidsID = x{1,1};

    %% 1. POWER (ELECTRODE SPACE)
    if ~exist(fullfile(params.power_folder,[bidsID '_power.mat']),'file')
        power = compute_power(params,bidsID);
        power.bidsID = bidsID;
        save(fullfile(params.power_folder,[bidsID '_power.mat']),'power')
    end
    % Plotting
    [power_fig, topoplot_fig] = plot_power(params,bidsID);
    saveas(power_fig,fullfile(params.figures_folder,[bidsID '_power.svg']));
    saveas(topoplot_fig,fullfile(params.figures_folder,[bidsID '_power_topoplots.svg']));
    close(power_fig);
    close(topoplot_fig);
    
    %% 2. PEAK FREQUENCY (ELECTRODE SPACE)
    if ~exist(fullfile(params.pf_folder,[bidsID '_peakfrequency.mat']),'file')
        pf = compute_peak_frequency(params,bidsID);
        pf.bidsID = bidsID;
        save(fullfile(params.pf_folder,[bidsID '_peakfrequency.mat']),'pf')
    end
    % Plotting
    pf_fig = plot_peakfrequency(params,bidsID);
    saveas(pf_fig,fullfile(params.figures_folder,[bidsID '_PeakFrequency.svg']));
    close(pf_fig)    
    
    %% Loop over frequency bands
    for iFreq = fields(params.freq_band)'
        
        %% 3. SOURCE RECONSTRUCTION (power at source space)
        if ~exist(fullfile(params.source_folder,[bidsID '_source_' iFreq{:} '.mat']),'file')
            source = compute_spatial_filter(params,bidsID,iFreq{:});
            source.bidsID = [bidsID '_' iFreq{:}];
            save(fullfile(params.source_folder,[bidsID '_source_' iFreq{:} '.mat']),'source')
        end
        
        %% 4.A FUNCTIONAL CONNECTIVITY - dwPLI
        if ~exist(fullfile(params.connectivity_folder,[bidsID '_dwpli_' iFreq{:} '.mat']),'file')
            connMatrix = compute_dwpli(params,bidsID,iFreq{:});
            save(fullfile(params.connectivity_folder,[bidsID '_dwpli_' iFreq{:} '.mat']),'connMatrix')
        end
        
        %% 4.B FUNCTIONAL CONNECTIVITY - AEC
        if ~exist(fullfile(params.connectivity_folder,[bidsID '_aec_' iFreq{:} '.mat']),'file')
            connMatrix = compute_aec(params,bidsID,iFreq{:});
            save(fullfile(params.connectivity_folder,[bidsID '_aec_' iFreq{:} '.mat']),'connMatrix')
        end
        
        %% 5. NETWORK CHARACTERIZATION (GRAPH MEASURES)
        if ~exist(fullfile(params.graph_folder,[bidsID '_graph_dwpli_' iFreq{:} '.mat']),'file')
            graph_dwpli = compute_graph_measures(params,bidsID,iFreq{:},'dwpli');
            save(fullfile(params.graph_folder,[bidsID '_graph_dwpli_' iFreq{:} '.mat']),'graph_dwpli')
        end        
        if ~exist(fullfile(params.graph_folder,[bidsID '_graph_aec_' iFreq{:} '.mat']),'file')
            graph_aec = compute_graph_measures(params,bidsID,iFreq{:},'aec');
            save(fullfile(params.graph_folder,[bidsID '_graph_aec_' iFreq{:} '.mat']),'graph_aec')
        end
        
    end
    % Plot connectivity matrices in all frequency bands and save it in figures folder
    dwpli_fig = plot_connectivity(params,bidsID,'dwpli');
    aec_fig = plot_connectivity(params,bidsID,'aec');
    saveas(dwpli_fig,fullfile(params.figures_folder,[bidsID '_dwpli.svg']));
    saveas(aec_fig,fullfile(params.figures_folder,[bidsID '_aec.svg']));
    close(dwpli_fig);
    close(aec_fig);

end
end

