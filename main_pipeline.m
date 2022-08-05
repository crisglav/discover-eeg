% Imports EEG dataset, preprocesses it, and extracts brain features
% 
% Cristina Gil, TUM, cristina.gil@tum.de, 25.07.2022
function main_pipeline()

clear all; close all;
rng('default'); % For reproducibility - See discussion in https://sccn.ucsd.edu/pipermail/eeglablist/2022/016932.html

% Define the parameters
params = define_params();
cd(params.main_folder)
save(fullfile(params.preprocessed_data_path,'pipeline_params.mat'),'params');

%% ======= IMPORT RAW DATA =========
% Import raw data in BIDS format
[STUDY, ALLEEG, bids] = pop_importbids(params.raw_data_path,'outputdir',params.preprocessed_data_path,...
    'studyName',params.study,'bidstask',params.task,'bidschanloc',params.bidschanloc,'bidsevent','off');

for iRec=1:length(ALLEEG)
    % Retrieve data
    EEGtemp = eeg_checkset(ALLEEG(iRec),'loaddata');
    
    % Add reference electrode
    EEGtemp = pop_chanedit(EEGtemp, 'append',EEGtemp.nbchan, ...
        'changefield', {EEGtemp.nbchan+1,'labels',bids.data(iRec).EEGReference},...
        'changefield', {EEGtemp.nbchan+1, 'X', params.RefCoord.X}, ...
        'changefield', {EEGtemp.nbchan+1, 'Y', params.RefCoord.Y}, ...
        'changefield', {EEGtemp.nbchan+1, 'Z', params.RefCoord.Z},...
        'setref',{['1:' num2str(EEGtemp.nbchan)],bids.data(iRec).EEGReference});
    
    % Use electrode positions from the electrodes.tsv file or from a standard template in the MNI coordinate system
    if strcmp(params.bidschanloc, 'on')
        % If electrode positions are chosen from the .tsv the coordinate
        % system might need to be adjusted (user has to define it in define_params.m)
        EEGtemp = pop_chanedit(EEGtemp, 'nosedir',params.nosedir);
        eegchans = find(contains(lower({ALLEEG(1).chanlocs.type}),'eeg'));
    else
        % Look for electrode positions in a standard template
        EEGtemp=pop_chanedit(EEGtemp, 'lookup','standard_1005.elc');
        non_standard_chans = cellfun(@isempty,{EEGtemp.chanlocs.X});
        eegchans = find(~non_standard_chans);
        if any(non_standard_chans)
            clabels = {EEGtemp.chanlocs(non_standard_chans).labels};
            c = sprintf('%s ', clabels{:});
            warning(['The position of the channel(s) ' c 'was not found in a standard template and they will be removed. If you want to include them please specify their position in a electrodes.tsv and change define_params accordingly.']);           
        end
    end
    
    % Select only EEG channels for preprocessing
    EEGtemp = pop_select(EEGtemp, 'channel', eegchans);
    EEGtemp.chaninfo.removedchans = [];
    
    % Save datafile, clear it from memory and store it in the ALLEEG structure
    EEGtemp = pop_saveset(EEGtemp, 'savemode', 'resave');
    EEGtemp.data = 'in set file';
    ALLEEG = eeg_store(ALLEEG, EEGtemp, iRec);
end  

% % OPTIONAL - Check that the electrodes positions are ok
% figure; topoplot([],ALLEEG(1).chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo',ALLEEG(1).chaninfo);
% hold on,
% figure; topoplot([],ALLEEG(1).chaninfo.nodatchans, 'style', 'blank',  'electrodes', 'labelpoint');

CURRENTSTUDY = 1;
EEG = ALLEEG;
CURRENTSET = 1:length(EEG);

%% ======== PREPROCESSING =========
tic
for iRec=1:length(EEG)
    % Retrieve data
    EEGtemp = eeg_checkset(EEG(iRec),'loaddata');
    
    % 1. CLEAN LINE NOISE
    EEGtemp = pop_cleanline(EEGtemp,'linefreqs',50,'newversion',1);
    
      
    % 2. REMOVE BAD CHANNELS
    [EEGtemp.urchanlocs] = deal(EEGtemp.chanlocs); % Keep original channels
    EEGtemp = pop_clean_rawdata(EEGtemp,'FlatlineCriterion', params.FlatlineCriterion,...
                                'ChannelCriterion',params.ChannelCriterion,...
                                'LineNoiseCriterion',params.LineNoiseCriterion,...
                                'Highpass',params.Highpass,...
                                'FuseChanRej',params.FuseChanRej,... 
                                'BurstCriterion','off',...
                                'WindowCriterion','off',...
                                'BurstRejection','off',...
                                'Distance','Euclidian',...
                                'WindowCriterionTolerances','off');
                            
    % 3. REREFERENCE TO AVERAGE REFERENCE
    if strcmp(params.addRefChannel,'on')
        EEGtemp = pop_reref(EEGtemp,[],'interpchan',[],'refloc', EEG(1).chaninfo.nodatchans);
    else
        EEGtemp = pop_reref(EEGtemp,[],'interpchan',[]);
    end

    % 4. REMOVE ARTIFACTS WITH ICA
    EEGtemp = pop_runica(EEGtemp,'icatype','runica','concatcond','off');
    EEGtemp = pop_iclabel(EEGtemp,'default');
    EEGtemp = pop_icflag(EEGtemp, params.IClabel); % flag artifactual components using IClabel
    classifications = EEGtemp.etc.ic_classification.ICLabel.classifications; % Keep classifications before component substraction
    EEGtemp = pop_subcomp(EEGtemp,[],0); % Substract artifactual independent components
    EEGtemp.etc.ic_classification.ICLabel.orig_classifications = classifications;


    % 5. INTERPOLATE MISSING CHANNELS
    urchanlocs = EEGtemp.urchanlocs;
    [~, iref] = setdiff({EEGtemp.chanlocs.labels},{EEGtemp.urchanlocs.labels}); % Handle reference in case it was added back
    if ~isempty(iref)
        urchanlocs(end+1) = EEGtemp.chanlocs(iref);
    end
    EEGtemp = pop_interp(EEGtemp, urchanlocs, 'spherical');
    

    % 6. REMOVE BAD TIME SEGMENTS
    EEGtemp = pop_clean_rawdata(EEGtemp,'FlatlineCriterion','off',...
                                'ChannelCriterion','off',...
                                'LineNoiseCriterion','off',...
                                'Highpass','off',...
                                'BurstCriterion',params.BurstCriterion,...
                                'WindowCriterion',params.WindowCriterion,...
                                'BurstRejection','on',...
                                'Distance','Euclidian',...
                                'WindowCriterionTolerances',params.WindowCriterionTolerances);
    EEGtemp.etc.eventsAfterCRD = EEGtemp.event; % Keep events for visualization later on.

    
    % 7. SEGMENT DATA INTO EPOCHS
    % EEGLab pop_epoch is designed to trim the data based on events. For
    % resting-state data, in which no events are defined, this function is
    % tricky to use. I added markers called 'epoch_start' each 10*(1-0.5) = 5 seconds with a duration of 1 sample.
    % Note: Epochs containing discontinutities will be automatically rejected.
    % Create markers each x seconds and add them at the end of existing event markers
    epoch_start = 1 : EEGtemp.srate * params.epoch_length * (1-params.epoch_overlap) : EEGtemp.pnts;
    nEpochs = length(epoch_start);
    if ~isempty(EEGtemp.event)
        events = EEGtemp.event;
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
    EEGtemp.event = events;
    EEGtemp = pop_epoch(EEGtemp,{'epoch_start'},[0 params.epoch_length],'epochinfo','yes');
    
    
    % Save datafile, clear it from memory and store it in the ALLEEG structure
    EEGtemp = pop_saveset(EEGtemp, 'savemode', 'resave');
    EEGtemp.data = 'in set file';
    EEG = eeg_store(EEG, EEGtemp, iRec);
end


% Save study
EEG = eeg_checkset(EEG);
EEG = pop_saveset(EEG,'savemode','resave');
ALLEEG = EEG;
STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', [params.study '_preprocessed'], 'filepath', params.preprocessed_data_path);

% PLOTTING PREPROCESSING
% Visualization of detected bad channels. If you set 'FuseChanRej' on, the union of bad channels in all tasks is rejected!
plot_badchannels(params,EEG);
% Visualization of rejected ICs
plot_ICs(params,EEG);
% Visualization of rejected time segments
plot_badtimesegments(params,EEG);

clear EEG ALLEEG;
params.preptime = toc;

%% ======= EXTRACTION OF BRAIN FEATURES =========

% You can start directly with preprocessed data in BIDS format by loading an EEGLAB STUDY
params = define_params();
cd(params.main_folder)
[STUDY, ~] = pop_loadstudy('filename', [params.study '_preprocessed.study'], 'filepath', params.preprocessed_data_path);

% % OPTIONAL - Visualization of corregistration of electroes and sources for one exemplary dataset (check that
% % electrodes are aligned with the head model)
% plot_electrodesandsources(params,'sub-01_task-closed')
% 
% % OPTIONAL -  Visualization of atlas regions by network
% plot_atlasregions(params);

tic
for iRec=1:length(STUDY.datasetinfo)
    
    % BIDS ID
    x = strsplit(STUDY.datasetinfo(iRec).filename,'_eeg.set');
    bidsID = x{1,1};

    % 1. POWER (ELECTRODE SPACE)
    if ~exist(fullfile(params.power_folder,[bidsID '_power.mat']),'file')
        compute_power(params,bidsID);
    end
    
    % 2. PEAK FREQUENCY (ELECTRODE SPACE)
    if ~exist(fullfile(params.power_folder,[bidsID '_peakfrequency.mat']),'file')
        compute_peakfrequency(params,bidsID);
    end    
    
    % Loop over frequency bands
    for iFreq = fields(params.freq_band)'
        
        % 3. SOURCE RECONSTRUCTION (power at source space)
        if ~exist(fullfile(params.source_folder,[bidsID '_source_' iFreq{:} '.mat']),'file')
            compute_spatial_filter(params,bidsID,iFreq{:});
        end
        
        % 4.A FUNCTIONAL CONNECTIVITY - dwPLI
        if ~exist(fullfile(params.connectivity_folder,[bidsID '_dwpli_' iFreq{:} '.mat']),'file')
            compute_dwpli(params,bidsID,iFreq{:});
        end
        
        % 4.B FUNCTIONAL CONNECTIVITY - AEC
        if ~exist(fullfile(params.connectivity_folder,[bidsID '_aec_' iFreq{:} '.mat']),'file')
            compute_aec(params,bidsID,iFreq{:});
        end
        
        % 5. NETWORK CHARACTERIZATION (GRAPH MEASURES)
        if ~exist(fullfile(params.graph_folder,[bidsID '_graph_dwpli_' iFreq{:} '.mat']),'file')
            compute_graph_measures(params,bidsID,iFreq{:},'dwpli');
        end        
        if ~exist(fullfile(params.graph_folder,[bidsID '_graph_aec_' iFreq{:} '.mat']),'file')
            compute_graph_measures(params,bidsID,iFreq{:},'aec');
        end
        
    end
    
    % PLOTTING
    % Plot source power in all frequency bands
    fig = plot_power_source(params,bidsID);
    saveas(fig,fullfile(params.source_folder,[bidsID '_source.svg']));
    close(fig);
    
    for iConMeas = {'dwpli','aec'}     
        % Plot connectivity matrices in all frequency bands and save them in connectivity folder
        fig = plot_connectivity(params,bidsID,iConMeas{:});
        saveas(fig,fullfile(params.connectivity_folder,[bidsID '_' iConMeas{:} '.svg']));
        close(fig);
        
        % Plot graph measures in all frequency bands and save them in the graph measures folder
        [f_degree, f_cc, f_global] = plot_graph_measures(params,bidsID,iConMeas{:});
        saveas(f_degree,fullfile(params.graph_folder,[bidsID '_' iConMeas{:} '_degree.svg']));
        saveas(f_cc,fullfile(params.graph_folder,[bidsID '_' iConMeas{:} '_cc.svg']));
        saveas(f_global,fullfile(params.graph_folder,[bidsID '_' iConMeas{:} '_global.svg']));
        close(f_degree, f_cc, f_global);   
    end
    
    % Generate individual recording reports with figures
    recording_report(params,bidsID)    
end
params.feattime = toc;
end



