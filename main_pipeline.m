% Imports EEG dataset, preprocesses it, and extracts brain features
% 
% Cristina Gil, TUM, cristina.gil@tum.de, 25.07.2022


clear all; close all;
rng('default'); % For reproducibility - See discussion in https://sccn.ucsd.edu/pipermail/eeglablist/2022/016932.html

% Define the parameters
params = define_params();
save(fullfile(params.preprocessed_data_path,'pipeline_params.mat'),'params');
 
if(isempty(gcp('nocreate')))
    parObj = parpool();
end

%% ======= IMPORT RAW DATA =========
% Try to load the already created study, otherwise import raw data with pop_importbids
if exist(fullfile(params.preprocessed_data_path,[params.study '.study']),'file')
    [STUDY, ALLEEG] = pop_loadstudy('filename', [params.study '.study'], 'filepath', params.preprocessed_data_path);
else
    % Import raw data in BIDS format
    [STUDY, ALLEEG] = pop_importbids(params.raw_data_path,'outputdir',params.preprocessed_data_path,...
        'studyName',params.study,'sessions',params.session,'runs',params.runs,'bidstask',params.task,...
        'bidschanloc',params.bidschanloc,'bidsevent','off');
    
    for iRec=1:length(ALLEEG)
        % Retrieve data
        EEGtemp = eeg_checkset(ALLEEG(iRec),'loaddata');
        
        % Add reference electrode
        EEGtemp = pop_chanedit(EEGtemp, 'append',EEGtemp.nbchan, ...
            'changefield', {EEGtemp.nbchan+1,'labels',ALLEEG(iRec).BIDS.tInfo.EEGReference},...
            'changefield', {EEGtemp.nbchan+1, 'X', params.RefCoord.X}, ...
            'changefield', {EEGtemp.nbchan+1, 'Y', params.RefCoord.Y}, ...
            'changefield', {EEGtemp.nbchan+1, 'Z', params.RefCoord.Z},...
            'setref',{['1:' num2str(EEGtemp.nbchan)],ALLEEG(iRec).BIDS.tInfo.EEGReference});
        
        % Use electrode positions from the electrodes.tsv file or from a standard template in the MNI coordinate system
        if strcmp(params.bidschanloc, 'on')
            % If electrode positions are chosen from the .tsv the coordinate
            % system might need to be adjusted (user has to define it in define_params.m)
            EEGtemp = pop_chanedit(EEGtemp, 'nosedir',params.nosedir);
            eegchans = find(contains(lower({ALLEEG(iRec).chanlocs.type}),'eeg'));
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
        STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', params.study, 'filepath', params.preprocessed_data_path);
        
    end
end
% % OPTIONAL - Check that the electrodes positions are ok
% figure; topoplot([],ALLEEG(1).chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo',ALLEEG(1).chaninfo);
% figure; topoplot([],ALLEEG(1).chaninfo.nodatchans, 'style', 'blank',  'electrodes', 'labelpoint');

%% ======== PREPROCESSING =========
% Find the latest preprocessed recording and start with the next one
first = find(cellfun(@isempty, {ALLEEG.setname}),1);
to_delete = {};

% Loop over the recordings that have not been preprocessed
for iRec=first:length(ALLEEG)

    % Retrieve data
    EEGtemp = eeg_checkset(ALLEEG(iRec),'loaddata');
    
    % OPTIONAL. DOWNSAMPLE DATA
    EEGtemp = pop_resample(EEGtemp, params.sampling_rate);
    
    % 1. CLEAN LINE NOISE
    try
        EEGtemp = pop_cleanline(EEGtemp,'linefreqs',EEGtemp.BIDS.tInfo.PowerLineFrequency,'newversion',1);
    catch ME
        warning(['CleanLine not performed: ' ME.message ' Make sure you specify the Line Noise Frequency in the *_eeg.json file'])
    end
    
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
    if(isfield(EEGtemp.etc,'clean_channel_mask'))
        clean_channel_mask = EEGtemp.etc.clean_channel_mask;
    end
    
    % 3. REREFERENCE TO AVERAGE REFERENCE
    if strcmp(params.addRefChannel,'on')
        EEGtemp = pop_reref(EEGtemp,[],'interpchan',[],'refloc', EEGtemp.chaninfo.nodatchans);
    else
        EEGtemp = pop_reref(EEGtemp,[],'interpchan',[]);
    end
    
    % ICA is non-deterministic and several runs of the algorithm lead to
    % small differences in the bad segment detection. We permorn 10 times
    % the ICA, interpolation of bad channels and bad segment detection and 
    % we select the run that is closer to the 'average' bad segment mask
    EEGOrig = EEGtemp;
    nRep = 10;
    EEGtemp_clean = cell(1,nRep);
    parfor iRep =1:nRep
        EEGtemp = EEGOrig;
        % 4. REMOVE ARTIFACTS WITH ICA
        EEGtemp = pop_runica(EEGtemp,'icatype','runica','concatcond','off');
        EEGtemp = pop_iclabel(EEGtemp,'default');
        EEGtemp = pop_icflag(EEGtemp, params.IClabel); % flag artifactual components using IClabel
        classifications = EEGtemp.etc.ic_classification.ICLabel.classifications; % Keep classifications before component substraction
        EEGtemp = pop_subcomp(EEGtemp,[],0); % Subtract artifactual independent components
        EEGtemp.etc.ic_classification.ICLabel.orig_classifications = classifications;
        
        % 5. INTERPOLATE MISSING CHANNELS
        urchanlocs = EEGtemp.urchanlocs;
        [~, iref] = setdiff({EEGtemp.chanlocs.labels},{EEGtemp.urchanlocs.labels}); % Handle reference in case it was added back
        if ~isempty(iref)
            urchanlocs(end+1) = EEGtemp.chanlocs(iref);
        end
        EEGtemp = pop_interp(EEGtemp, urchanlocs, 'spherical');
        
        
        % 6. REMOVE BAD TIME SEGMENTS
        EEGtemp_dirty{iRep} = EEGtemp;
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
        EEGtemp_clean{iRep} = EEGtemp;
    end
        
    % Select the run closer to the 'average' bad time segments mask
    etc = cellfun(@(x) x.etc, EEGtemp_clean, 'UniformOutput',0);
    csm = cell2mat(cellfun(@(x) x.clean_sample_mask, etc, 'UniformOutput',0)');
    csmAverage = mean(csm,1);
    distToAvg = sum(abs(csm - csmAverage), 2);
    [~,selRun] = min(distToAvg);
    EEGtemp = EEGtemp_clean{selRun};
    clear EEGtemp_clean
    
    % 7. SEGMENT DATA INTO EPOCHS
    % EEGLab pop_epoch is designed to trim the data based on events. For
    % resting-state data, in which no events are defined, this function is
    % tricky to use. I added markers called 'epoch_start' each 2*(1-0.5) = 1 second with a duration of 1 sample.
    % Note: Epochs containing discontinutities will be automatically rejected.
    % Create markers each x seconds and add them at the end of existing event markers
    try
        % Create spatially distributed markers
        EEGtemp = eeg_regepochs(EEGtemp,'recurrence',params.epoch_length * (1-params.epoch_overlap),'eventtype','epoch_start','extractepochs','off');
        % Segment the data into epochs
        EEGtemp = pop_epoch(EEGtemp,{'epoch_start'},[0 params.epoch_length],'epochinfo','yes'); % epoch latencies are lost in this step
        % Mark recordings without any trial left to remove later on
        if EEGtemp.trials == 0, to_delete{end +1} = EEGtemp.filename; end
    catch ME
        warning(['Data segmentation not performed: ' ME.message ' Probably no clean epochs remained. Recording will be removed.'])
        to_delete{end +1} = EEGtemp.filename;
    end
    
    
    % Save datafile, clear it from memory and store it in the ALLEEG structure
    if exist('clean_channel_mask','var'), EEGtemp.etc.clean_channel_mask = clean_channel_mask; end; clear 'clean_channel_mask';
    EEGtemp = pop_saveset(EEGtemp, 'savemode', 'resave');
    EEGtemp.data = 'in set file';
    ALLEEG = eeg_store(ALLEEG, EEGtemp, iRec);
    STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', params.study, 'filepath', params.preprocessed_data_path);

end

% PLOTTING PREPROCESSING
% Visualization of detected bad channels. If you set 'FuseChanRej' on, the union of bad channels in all tasks is rejected!
plot_badchannels(params,ALLEEG);
% Visualization of rejected ICs
plot_ICs(params,ALLEEG);
% Visualization of rejected time segments
plot_badtimesegments(params,ALLEEG);
% Visualize number of clean epochs per recording
plot_epochs(params,ALLEEG);

% Delete recordings with no epochs remaining
mask = matches({ALLEEG.filename},to_delete);
ALLEEG(mask) = [];
STUDY.datasetinfo(mask) = [];
s = split(to_delete,{'_'});
s = unique(s(:,:,1));
mask = matches(STUDY.subject,s);
STUDY.subject(mask) = [];

% Save study
STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', [params.study '-clean'], 'filepath', params.preprocessed_data_path);

clear EEGtemp ALLEEG;

%% ======= EXTRACTION OF BRAIN FEATURES =========

% % You can start directly with preprocessed data in BIDS format by loading an EEGLAB STUDY
% params = define_params();
% [STUDY, ~] = pop_loadstudy('filename', [params.study '-clean.study'], 'filepath', params.preprocessed_data_path);
%%
% % OPTIONAL - Visualization of corregistration of electroes and sources for one exemplary dataset (check that
% % electrodes are aligned with the head model)
% plot_electrodesandsources(params,'sub-001_task-EC')
% 
% % OPTIONAL -  Visualization of atlas regions by network
% plot_atlasregions(params);
freqs = fields(params.freq_band);
conMeas = {'dwpli','aec'};
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
    for iFreq = 1:length(freqs)
        
        % 3. SOURCE RECONSTRUCTION (power at source space)
        if ~exist(fullfile(params.source_folder,[bidsID '_source_' freqs{iFreq} '.mat']),'file')
            compute_spatial_filter(params,bidsID,freqs{iFreq});
        end
        
        % 4.A FUNCTIONAL CONNECTIVITY - dwPLI
        if ~exist(fullfile(params.connectivity_folder,[bidsID '_dwpli_' freqs{iFreq} '.mat']),'file')
            try
                compute_dwpli(params,bidsID,freqs{iFreq});
            catch ME
                warning([bidsID ' - ' ME.message]);
                continue;
            end
        end
       
        % 4.B FUNCTIONAL CONNECTIVITY - AEC
        if ~exist(fullfile(params.connectivity_folder,[bidsID '_aec_' freqs{iFreq} '.mat']),'file')
            try
                compute_aec(params,bidsID,freqs{iFreq});
            catch ME
                warning([bidsID ' - ' ME.message]);
                continue;
            end
        end
        
        % 5. NETWORK CHARACTERIZATION (GRAPH MEASURES)
        if ~exist(fullfile(params.graph_folder,[bidsID '_graph_dwpli_' freqs{iFreq} '.mat']),'file')
            try
                compute_graph_measures(params,bidsID,freqs{iFreq},'dwpli');
            catch ME
                warning([bidsID ' - ' ME.message]);
                continue;
            end
        end
        if ~exist(fullfile(params.graph_folder,[bidsID '_graph_aec_' freqs{iFreq} '.mat']),'file')
            try
                compute_graph_measures(params,bidsID,freqs{iFreq},'aec');
            catch ME
                warning([bidsID ' - ' ME.message]);
                continue;
            end
        end
        
    end
    
    % PLOTTING
    % Plot source power in all frequency bands
    if ~exist(fullfile(params.source_folder,[bidsID '_source.svg']),'file')
        fig = plot_power_source(params,bidsID);
        saveas(fig,fullfile(params.source_folder,[bidsID '_source.svg']));
        close(fig);
    end
    
    for iConMeas = 1:length(conMeas)
        % Plot connectivity matrices in all frequency bands and save them in connectivity folder
        if ~exist(fullfile(params.connectivity_folder,[bidsID '_' conMeas{iConMeas} '.svg']),'file')
            try
                fig = plot_connectivity(params,bidsID,conMeas{iConMeas});
                saveas(fig,fullfile(params.connectivity_folder,[bidsID '_' conMeas{iConMeas} '.svg']));
                close(fig);
            catch ME
                warning([bidsID ' - ' ME.message]);
            end
        end

        
        % Plot graph measures in all frequency bands and save them in the graph measures folder
        if ~exist(fullfile(params.graph_folder,[bidsID '_' conMeas{iConMeas} '_degree.svg']),'file') ||...
                ~exist(fullfile(params.graph_folder,[bidsID '_' conMeas{iConMeas} '_cc.svg']),'file') ||...
                ~exist(fullfile(params.graph_folder,[bidsID '_' conMeas{iConMeas} '_global.svg']),'file')
            try
                [f_degree, f_cc, f_global] = plot_graph_measures(params,bidsID,conMeas{iConMeas});
                saveas(f_degree,fullfile(params.graph_folder,[bidsID '_' conMeas{iConMeas} '_degree.svg']));
                saveas(f_cc,fullfile(params.graph_folder,[bidsID '_' conMeas{iConMeas} '_cc.svg']));
                saveas(f_global,fullfile(params.graph_folder,[bidsID '_' conMeas{iConMeas} '_global.svg']));
                close(f_degree, f_cc, f_global);
            catch ME
                warning([bidsID ' - ' ME.message]);
            end
                
        end
    end
    
    % Generate individual recording reports with figures
    if ~exist(fullfile(params.reports_folder,[bidsID '_report.pdf']),'file')
        try
            recording_report(params,bidsID);
        catch ME
            warning([bidsID ' - ' ME.message]);
        end
        
    end
end




