% Imports EEG dataset, preprocesses it, and extracts brain features
% 
% Cristina Gil, TUM, cristina.gil@tum.de, 25.07.2022

clear all; close all;
rng('default'); % For reproducibility - See discussion in https://sccn.ucsd.edu/pipermail/eeglablist/2022/016932.html

% Define the parameters
params = define_params('params.json');

%% ======= IMPORT RAW DATA =========
% Try to load the already created study, otherwise import raw data with pop_importbids
if exist(fullfile(params.PreprocessedDataPath,[params.StudyName '.study']),'file')
    [STUDY, ALLEEG] = pop_loadstudy('filename', [params.StudyName '.study'], 'filepath', params.PreprocessedDataPath);
else
    % Import raw data in BIDS format
    [STUDY, ALLEEG] = pop_importbids(params.RawDataPath,'outputdir',params.PreprocessedDataPath,...
        'studyName',params.StudyName,'sessions',params.Session,'runs',params.Run,'bidstask',params.Task,...
        'bidschanloc',params.BidsChanloc,'bidsevent','off');
    
    for iRec=1:length(ALLEEG)
        % Retrieve data
        EEGtemp = eeg_checkset(ALLEEG(iRec),'loaddata');
        
        % For event-related data, check if the recording has the specified EventMarker
        if params.PreprocEventData
            if any(strcmp({EEGtemp.event.type},params.EventMarker))
                fprintf('We detected %d instances of marker %s.\n',sum(strcmp({EEGtemp.event.type},params.EventMarker)), params.EventMarker);
            else
                error('We could not find EventMarker in the data.')
            end
        end
        
        % Add reference electrode
        EEGtemp = pop_chanedit(EEGtemp, 'append',EEGtemp.nbchan, ...
            'changefield', {EEGtemp.nbchan+1,'labels',ALLEEG(iRec).BIDS.tInfo.EEGReference},...
            'changefield', {EEGtemp.nbchan+1, 'X', params.RefCoord.X}, ...
            'changefield', {EEGtemp.nbchan+1, 'Y', params.RefCoord.Y}, ...
            'changefield', {EEGtemp.nbchan+1, 'Z', params.RefCoord.Z},...
            'setref',{['1:' num2str(EEGtemp.nbchan)],ALLEEG(iRec).BIDS.tInfo.EEGReference});
        
        % Use electrode positions from the electrodes.tsv file or from a standard template in the MNI coordinate system
        if strcmp(params.BidsChanloc, 'on')
            % If electrode positions are chosen from the .tsv, the coordinate
            % system might need to be adjusted (user has to define it in params.json)
            EEGtemp = pop_chanedit(EEGtemp, 'nosedir',params.NoseDir);
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
        STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', params.StudyName, 'filepath', params.PreprocessedDataPath);
        
    end
end
% % OPTIONAL - Check that the electrode positions are ok
% figure; topoplot([],ALLEEG(1).chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo',ALLEEG(1).chaninfo);
% figure; topoplot([],ALLEEG(1).chaninfo.nodatchans, 'style', 'blank',  'electrodes', 'labelpoint','chaninfo',ALLEEG(1).chaninfo);

%% ======== PREPROCESSING =========
% Find the latest preprocessed recording and start with the next one
not_preprocessed = find(~cellfun(@(x) strcmp(x,'preprocessed'), {ALLEEG.comments}));
to_delete = {};

% Log file of all the recordings
fid0 = fopen(fullfile(params.ReportsPath,'preprocessing.log'),'a');

% Loop over the recordings that have not been preprocessed
for iRec = not_preprocessed
    
    t1 = tic;
    % Log file of single recordings
    id = regexp(ALLEEG(iRec).filename,'.*(?=(_eeg.set))','match');
    id = id{:};
    fid = fopen(fullfile(params.ReportsPath,[id '.log']),'a');

    % Retrieve data
    EEGtemp = eeg_checkset(ALLEEG(iRec),'loaddata');
    
    % OPTIONAL. DOWNSAMPLE DATA
    if(~isempty(params.DownsamplingRate))
    try
        EEGtemp = pop_resample(EEGtemp, params.DownsamplingRate);
        fprintf(fid,'0. Downsampling performed. \n');
    catch ME
        fprintf(fid,['--- Downsampling not performed: ' ME.message '. \n']);
    end
    end
    
    % 1. CLEAN LINE NOISE
    try
        EEGtemp = pop_cleanline(EEGtemp,'linefreqs',EEGtemp.BIDS.tInfo.PowerLineFrequency,'newversion',1);
        fprintf(fid,'1. CleanLine performed. \n');
    catch ME
        fprintf(fid,['--- CleanLine not performed: ' ME.message ' Make sure you specify the Line Noise Frequency in the *_eeg.json file. \n']);
        fprintf(fid0,'--- %s not preprocessed successfully.\n', id);
    end
    
    % 2. REMOVE BAD CHANNELS
    try
        [EEGtemp.urchanlocs] = deal(EEGtemp.chanlocs); % Keep original channels
        EEGtemp = pop_clean_rawdata(EEGtemp,'FlatlineCriterion', params.FlatLineCriterion,...
            'ChannelCriterion',params.ChannelCriterion,...
            'LineNoiseCriterion',params.LineNoiseCriterion,...
            'Highpass',params.HighPass,...
            'BurstCriterion','off',...
            'WindowCriterion','off',...
            'BurstRejection','off',...
            'Distance','Euclidian',...
            'WindowCriterionTolerances','off');
        if(isfield(EEGtemp.etc,'clean_channel_mask'))
            clean_channel_mask = EEGtemp.etc.clean_channel_mask;
        end
        fprintf(fid,'2. Bad channel removal performed. \n');
    catch ME
        fprintf(fid, ['--- Bad channel removal not performed: ' ME.message '. \n']);
        fprintf(fid0,'--- %s not preprocessed successfully.\n', id);
    end
    
    % 3. REREFERENCE TO AVERAGE REFERENCE
    % Note: bad channels from point 2 are interpolated before re-refencing to
    % the average, but these interpolated channels are not carried along.
    try
        if strcmp(params.AddRefChannel,'on')
            EEGtemp = pop_reref(EEGtemp,[],'interpchan',[],'refloc', EEGtemp.chaninfo.nodatchans);
        else
            EEGtemp = pop_reref(EEGtemp,[],'interpchan',[]);
        end
        % Do not track interpolated electrodes for re-referencing
        EEGtemp.chaninfo.removedchans = EEGtemp.chaninfo.removedchans(~ismember({EEGtemp.chaninfo.removedchans.ref}, 'average'));
        fprintf(fid,'3. Re-referencing performed. \n');
    catch ME
        fprintf(fid,['--- Re-referencing not performed: ' ME.message '. \n']);
        fprintf(fid0,'--- %s not preprocessed successfully.\n', id);
    end
    
    % SPECIFIC FOR EVENT-RELATED DATA - if you are only interested in
    % preprocessing certain segments surrounding an event, but not in time intervals in
    % between those segments (e.g. because they contain artifacts), you can
    % cut out these segments of interest here. They are concatenated and 
    % then fed into the ICA. 
    if params.PreprocEventData
        EEGtemp = pop_epoch(EEGtemp,{params.EventMarker},params.EventBounds);
        EEGtemp = eeg_epoch2continuous(EEGtemp);
    end    
    
    % 4. REMOVE ARTIFACTS WITH ICA
    % 5. INTERPOLATE MISSING CHANNELS
    % 6. REMOVE BAD TIME SEGMENTS

    % ICA is non-deterministic and several runs of the algorithm lead to
    % small differences in the bad segment detection. We perform 10 times
    % the ICA, interpolation of bad channels and bad segment detection and 
    % we select the run that is closest to the 'average' bad segment mask
    error_parallelization = false;
    try
        nRep = params.NICARepetitions;
        EEGtemp_clean = cell(1,nRep);
        parfor iRep =1:nRep
            % Steps 4, 5, and 6
            EEGtemp_clean{iRep} = preprocessing_ICA(EEGtemp,params);
        end
               
        % Log warnings if some repetitions were not performed
        error_flags = cellfun(@isnumeric, EEGtemp_clean);
        error_flags = cell2mat(EEGtemp_clean(error_flags));
        if ~isempty(error_flags)
            fprintf(fid,'--- Warning: %d repetitions failed at ICA. \n', sum(error_flags == 4));
            fprintf(fid,'--- Warning: %d repetitions failed at channel interpolation. \n', sum(error_flags == 5));
            fprintf(fid,'--- Warning: %d repetitions failed at bad channel rejection. \n', sum(error_flags == 6));
        end
        
        % Custom function that selects the repetition closest to the
        % 'average' bad time segments mask
        EEGtemp = preprocessing_select_ICA_rep(EEGtemp_clean);
        
        % Log success
        fprintf(fid,'4. ICA performed. \n5. Channel interpolation performed. \n6. Bad interval removal performed. \n');                
    catch
        fprintf(fid,'--- No ICA repetition was successful with parallelization. Trying without parallelization... \n');
        error_parallelization = true;
    end
    
    % Try steps 4, 5, 6 without parallelization
    if (error_parallelization)
        try
            EEGtemp_clean = cell(1,nRep);
            for iRep =1:nRep
                % Steps 4, 5 and 6
                EEGtemp_clean{iRep} = preprocessing_ICA(EEGtemp,params);
            end
            
            % Log warnings if some repetitions were not performed
            error_flags = cellfun(@isnumeric, EEGtemp_clean);
            error_flags = cell2mat(EEGtemp_clean(error_flags));
            if ~isempty(error_flags)
                fprintf(fid,'--- Warning: %d repetitions failed at ICA. \n', sum(error_flags == 4));
                fprintf(fid,'--- Warning: %d repetitions failed at channel interpolation. \n', sum(error_flags == 5));
                fprintf(fid,'--- Warning: %d repetitions failed at bad channel rejection. \n', sum(error_flags == 6));
                fprintf(fid0,'--- %s Warning. Not all ICA repetitions were successful.\n', id);
            end
            
            % Custom function that selects the repetition closest to the
            % 'average' bad time segments mask
            EEGtemp = preprocessing_select_ICA_rep(EEGtemp_clean);
            
            % Log success
            fprintf(fid,'4. ICA performed. \n5. Channel interpolation performed. \n6. Bad interval removal performed. \n');        

        catch
            fprintf(fid,['--- No ICA repetition was successful without parallelization. \n Exclude recording.\n']);
            to_delete{end +1} = EEGtemp.filename;
            fprintf(fid0,'--- %s not preprocessed successfully.\n', id);
            continue; % to next recording
        end
    end
    
        
    % 7. SEGMENT DATA INTO EPOCHS
    % EEGLab pop_epoch is designed to trim the data based on events. For
    % resting-state data, in which no events are defined, this function is
    % tricky to use. I added markers called 'epoch_start' each 2*(1-0.5) = 1 second with a duration of 1 sample.
    % Note: Epochs containing discontinuities will be automatically rejected.
    % Create markers each x seconds and add them at the end of existing event markers
    if(~isempty(params.EpochLength))
    try
        % Create spatially distributed markers
        EEGtemp = eeg_regepochs(EEGtemp,'recurrence',params.EpochLength * (1-params.EpochOverlap),'eventtype','epoch_start','extractepochs','off');
        % Segment the data into epochs
        EEGtemp = pop_epoch(EEGtemp,{'epoch_start'},[0 params.EpochLength],'epochinfo','yes'); % epoch latencies are lost in this step
        % Mark recordings without any trial left to remove later on
        if EEGtemp.trials == 0, to_delete{end +1} = EEGtemp.filename; end
        fprintf(fid,'7. Segmentation into epochs performed. \n');
    catch ME
        fprintf(fid,['Segmentation into epochs not performed: ' ME.message ' Probably no clean epochs remained. \n Exclude recording.\n']);
        to_delete{end +1} = EEGtemp.filename;
        fprintf(fid0,'--- %s not preprocessed successfully.\n', id);
        continue;
    end
    end

    % Save datafile, clear it from memory, and store it in the ALLEEG structure
    EEGtemp.comments = 'preprocessed';
    if exist('clean_channel_mask','var'), EEGtemp.etc.clean_channel_mask = clean_channel_mask; end; clear 'clean_channel_mask';
    EEGtemp = pop_saveset(EEGtemp, 'savemode', 'resave');
    EEGtemp.data = 'in set file';
    ALLEEG = eeg_store(ALLEEG, EEGtemp, iRec);
    STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', params.StudyName, 'filepath', params.PreprocessedDataPath);
    
    fprintf(fid,'Preprocessed data saved to disk');
    fclose(fid);
    t2 = toc(t1);
    fprintf(fid0,'%s preprocessed successfully. It took %.2f seconds.\n', id, t2);

end
fclose(fid0);
% Delete recordings with no epochs remaining
mask = matches({ALLEEG.filename},to_delete);
ALLEEG(mask) = [];
STUDY.datasetinfo(mask) = [];
s = split(to_delete,{'_'});
s = unique(s(:,:,1));
mask = matches(STUDY.subject,s);
STUDY.subject(mask) = [];

% Save study
STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', [params.StudyName '-clean'], 'filepath', params.PreprocessedDataPath);

% PLOTTING PREPROCESSING
% Visualization of detected bad channels. If you set 'FuseChanRej' on, the union of bad channels in all tasks is rejected!
plot_badchannels(params,ALLEEG);
% Visualization of rejected ICs
plot_ICs(params,ALLEEG);
% Visualization of rejected time segments
plot_badtimesegments(params,ALLEEG);
% Visualize number of clean epochs per recording
if ~isempty(params.EpochLength)
    plot_epochs(params,ALLEEG);
end
clear EEGtemp ALLEEG;

%% ======= EXTRACTION OF BRAIN FEATURES =========

% % You can start directly with preprocessed data in BIDS format by loading an EEGLAB STUDY
% params = define_params('demo_data/derivatives_v2023_10_05/params.json');
% [STUDY, ALLEEG] = pop_loadstudy('filename', [params.StudyName '-clean.study'], 'filepath', params.PreprocessedDataPath);
%%
% % OPTIONAL - Visualization of corregistration of electrodes and sources for one exemplary dataset (to check that
% % electrodes are aligned with the head model)
% plot_electrodesandsources(params,'sub-010002')
% 
% % OPTIONAL -  Visualization of atlas regions by network
% plot_atlasregions(params);
if params.BrainFeatExtr
    
    freqs = fields(params.FreqBand);
    conMeas = {'dwpli','aec'};
    
    for iRec=1:length(STUDY.datasetinfo)
        
        % BIDS ID
        x = strsplit(STUDY.datasetinfo(iRec).filename,'_eeg.set');
        bidsID = x{1,1};
        
        % 1. POWER (ELECTRODE SPACE)
        if ~exist(fullfile(params.PowerPath,[bidsID '_power.mat']),'file')
            compute_power(params,bidsID);
        end
        
        % 2. PEAK FREQUENCY (ELECTRODE SPACE)
        if ~exist(fullfile(params.PowerPath,[bidsID '_peakfrequency.mat']),'file')
            compute_peakfrequency(params,bidsID);
        end
        
        % Loop over frequency bands
        for iFreq = 1:length(freqs)
            
            % 3. SOURCE RECONSTRUCTION (power at source space)
            if ~exist(fullfile(params.SourcePath,[bidsID '_source_' freqs{iFreq} '.mat']),'file')
                compute_spatial_filter(params,bidsID,freqs{iFreq});
            end
            
            % 4.A FUNCTIONAL CONNECTIVITY - dwPLI
            if ~exist(fullfile(params.ConnectivityPath,[bidsID '_dwpli_' freqs{iFreq} '.mat']),'file')
                try
                    compute_dwpli(params,bidsID,freqs{iFreq});
                catch ME
                    warning([bidsID ' - ' ME.message]);
                    continue;
                end
            end
            
            % 4.B FUNCTIONAL CONNECTIVITY - AEC
            if ~exist(fullfile(params.ConnectivityPath,[bidsID '_aec_' freqs{iFreq} '.mat']),'file')
                try
                    compute_aec(params,bidsID,freqs{iFreq});
                catch ME
                    warning([bidsID ' - ' ME.message]);
                    continue;
                end
            end
            
            % 5. NETWORK CHARACTERIZATION (GRAPH MEASURES)
            if ~exist(fullfile(params.GraphPath,[bidsID '_graph_dwpli_' freqs{iFreq} '.mat']),'file')
                try
                    compute_graph_measures(params,bidsID,freqs{iFreq},'dwpli');
                catch ME
                    warning([bidsID ' - ' ME.message]);
                    continue;
                end
            end
            if ~exist(fullfile(params.GraphPath,[bidsID '_graph_aec_' freqs{iFreq} '.mat']),'file')
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
        if ~exist(fullfile(params.SourcePath,[bidsID '_source.svg']),'file')
            fig = plot_power_source(params,bidsID);
            saveas(fig,fullfile(params.SourcePath,[bidsID '_source.svg']));
            close(fig);
        end
        
        for iConMeas = 1:length(conMeas)
            % Plot connectivity matrices in all frequency bands and save them in connectivity folder
            if ~exist(fullfile(params.ConnectivityPath,[bidsID '_' conMeas{iConMeas} '.svg']),'file')
                try
                    fig = plot_connectivity(params,bidsID,conMeas{iConMeas});
                    saveas(fig,fullfile(params.ConnectivityPath,[bidsID '_' conMeas{iConMeas} '.svg']));
                    close(fig);
                catch ME
                    warning([bidsID ' - ' ME.message]);
                end
            end
            
            
            % Plot graph measures in all frequency bands and save them in the graph measures folder
            if ~exist(fullfile(params.GraphPath,[bidsID '_' conMeas{iConMeas} '_degree.svg']),'file') ||...
                    ~exist(fullfile(params.GraphPath,[bidsID '_' conMeas{iConMeas} '_cc.svg']),'file') ||...
                    ~exist(fullfile(params.GraphPath,[bidsID '_' conMeas{iConMeas} '_global.svg']),'file')
                try
                    [f_degree, f_cc, f_global] = plot_graph_measures(params,bidsID,conMeas{iConMeas});
                    saveas(f_degree,fullfile(params.GraphPath,[bidsID '_' conMeas{iConMeas} '_degree.svg']));
                    close(f_degree);
                    saveas(f_cc,fullfile(params.GraphPath,[bidsID '_' conMeas{iConMeas} '_cc.svg']));
                    close(f_cc);
                    saveas(f_global,fullfile(params.GraphPath,[bidsID '_' conMeas{iConMeas} '_global.svg']));
                    close(f_global);
                catch ME
                    warning([bidsID ' - ' ME.message]);
                end
                
            end
        end
        
        %     % Generate individual recording reports with figures
        %     if ~exist(fullfile(params.ReportsPath,[bidsID '_report.pdf']),'file')
        %         try
        %             recording_report(params,bidsID);
        %         catch ME
        %             warning([bidsID ' - ' ME.message]);
        %         end
        %
        %     end
        
        delete(fullfile(params.PreprocessedDataPath,bidsID,'eeg',[bidsID '_eeg.mat']));
    end
end



