% Imports EEG dataset, preprocesses it, and extracts brain features
%
% Cristina Gil, TUM, cristina.gil@tum.de, 25.07.2022
function main_pipelineV2()

clear all; %close all;
rng('default'); % For reproducibility - See discussion in https://sccn.ucsd.edu/pipermail/eeglablist/2022/016932.html

% Define the parameters
params = define_params();
cd(params.main_folder)
save(fullfile(params.preprocessed_data_path,'pipeline_params.mat'),'params');

if(isempty(gcp('nocreate')))
    parObj = parpool(4);
end

%% ======= IMPORT RAW DATA =========
% Import raw data in BIDS format
[STUDY, ALLEEG] = pop_importbids(params.raw_data_path,'outputdir',params.preprocessed_data_path,...
    'studyName',params.study,'bidstask',params.task,'bidschanloc',params.bidschanloc,'bidsevent','off');


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
    STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', params.study, 'filepath', params.preprocessed_data_path);
    
end

% % OPTIONAL - Check that the electrodes positions are ok
% figure; topoplot([],ALLEEG(1).chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo',ALLEEG(1).chaninfo);
% figure; topoplot([],ALLEEG(1).chaninfo.nodatchans, 'style', 'blank',  'electrodes', 'labelpoint');

% parpool('local');
%% ======== PREPROCESSING =========
tic
to_delete = {};
for iRec=1:length(ALLEEG)
    
    % Retrieve data
    EEGtemp = eeg_checkset(ALLEEG(iRec),'loaddata');
    
    % OPTIONAL. DOWNSAMPLE DATA
    rmpath(fullfile(params.fieldtrip_path, 'external', 'signal'))
    EEGtemp = pop_resample(EEGtemp, params.sampling_rate);
    addpath(fullfile(params.fieldtrip_path, 'external', 'signal'))
    
    % 1. CLEAN LINE NOISE
    try
        EEGtemp = pop_cleanline(EEGtemp,'linefreqs',EEGtemp.BIDS.tInfo.PowerLineFrequency,'newversion',1);
    catch ME
        warning(['CleanLine not performed: ' ME.message ' Make sure you specify the Line Noise Frequency in the *_eeg.json file'])
    end
    
    % 2. REMOVE BAD CHANNELS
    [EEGtemp.urchanlocs] = deal(EEGtemp.chanlocs); % Keep original channels
    EEGtempX = pop_clean_rawdata(EEGtemp,'FlatlineCriterion', params.FlatlineCriterion,...
        'ChannelCriterion',params.ChannelCriterion,... % could aslo be set to 'off'
        'LineNoiseCriterion',params.LineNoiseCriterion,...
        'Highpass',params.Highpass,...
        'FuseChanRej',params.FuseChanRej,...
        'BurstCriterion','off',...
        'WindowCriterion','off',...
        'BurstRejection','off',...
        'Distance','Euclidian',...
        'WindowCriterionTolerances','off');
    
    if(isfield(EEGtempX.etc,'clean_channel_mask'))
        clean_channel_mask = EEGtempX.etc.clean_channel_mask;
    end
    EEGtemp = EEGtempX;
    
    % 3. REREFERENCE TO AVERAGE REFERENCE
    if strcmp(params.addRefChannel,'on')
        EEGtemp = pop_reref(EEGtemp,[],'interpchan',[],'refloc', ALLEEG(iRec).chaninfo.nodatchans);
    else
        EEGtemp = pop_reref(EEGtemp,[],'interpchan',[]);
    end
    
    
    nRepet = 10; % outer loop repetitions to demonstrate enhanced reproducibility
    paramsIClabel = params.IClabel;
    parfor iRepet = 1:nRepet
        
        % 4. REMOVE ARTIFACTS WITH ICA
        nRepetICA = 10; % inner loop repetitions to identify most representative outcome
        EEGtempX2 = cell(nRepetICA,1)
        EEGtempX3 = cell(nRepetICA,1)
        for iRepetICA = 1:nRepetICA 
            EEGtempZ = pop_runica(EEGtemp,'icatype','runica','concatcond','off', 'stop', 1e-7,'maxsteps',100);
            EEGtempY = pop_iclabel(EEGtempZ,'default');
            EEGtempY = pop_icflag(EEGtempY, paramsIClabel); % flag artifactual components using IClabel
            
            classifications = EEGtempY.etc.ic_classification.ICLabel.classifications; % Keep classifications before component substraction
            EEGtempY = pop_subcomp(EEGtempY,[],0); % Subtract artifactual independent components
            EEGtempY.etc.ic_classification.ICLabel.orig_classifications = classifications;
        
            % 5. INTERPOLATE MISSING CHANNELS
            urchanlocs = EEGtempY.urchanlocs;
            [~, iref] = setdiff({EEGtempY.chanlocs.labels},{EEGtempY.urchanlocs.labels}); % Handle reference in case it was added back
            if ~isempty(iref)
                urchanlocs(end+1) = EEGtempY.chanlocs(iref);
            end
            EEGtempY = pop_interp(EEGtempY, urchanlocs, 'spherical');
            EEGtempX2{iRepetICA} = EEGtempY;

            % 6. DETECT BAD TIME SEGMENTS
            % EEGtempOld = EEGtempX;  
            EEGtempX3{iRepetICA} = pop_clean_rawdata(EEGtempY,'FlatlineCriterion','off',...
                'ChannelCriterion','off',...
                'LineNoiseCriterion','off',...
                'Highpass','off',...
                'BurstCriterion',params.BurstCriterion,...
                'WindowCriterion',params.WindowCriterion,...
                'BurstRejection','on',...
                'Distance','Euclidian',...
                'WindowCriterionTolerances',params.WindowCriterionTolerances);
            
            % get bad segment mask
            if(iRepetICA == 1)
                csm = EEGtempX3{iRepetICA}.etc.clean_sample_mask;
                csmAll = zeros(nRepetICA, length(csm)); 
                csmAll(iRepetICA,:) = ~csm;
            else
                csmAll(iRepetICA,:) = ~EEGtempX3{iRepetICA}.etc.clean_sample_mask;
            end  
        end
        
        %csmIntersect = sum(csmAll) == size(csmAll,1);
        csmAverage = mean(csmAll,1);
        
        % identify ICA run that is closest to average result
        distToAvg = sum(abs(csmAll - csmAverage), 2);
        [~,selRun] = min(distToAvg);
        
        EEGselected{iRepet} = EEGtempX3{selRun};
        EEGselectedBefore{iRepet} = EEGtempX2{selRun};  
    end

    % visualize bad segments of various repetitions
    chanSet = [3];    
    for selChan = chanSet
        figure
        hold on;
        datMat = EEGselectedBefore{1}.data(selChan,:);
        plot(datMat)
        display(selChan)
        for iRepet = 1:nRepet
            segMask = EEGselected{iRepet}.etc.clean_sample_mask;
            blub = 1:length(segMask);
            datMat = EEGselectedBefore{iRepet}.data(selChan,:);
            plot(blub(~segMask), datMat(~segMask)+40*iRepet)
        end
    end
    
    
    for iRepet = 1:nRepet
        EEGtemp = EEGselected{iRepet};
        EEGtemp.etc.eventsAfterCRD = EEGtemp.event; % Keep events for visualization later on.
        EEGtemp0{iRepet} = EEGtemp;
        
        % 7. SEGMENT DATA INTO EPOCHS
        % EEGLab pop_epoch is designed to trim the data based on events. For
        % resting-state data, in which no events are defined, this function is
        % tricky to use. I added markers called 'epoch_start' each 10*(1-0.5) = 5 seconds with a duration of 1 sample.
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
            warning(['Data segmentation not performed: ' ME.message ' Probably no clean epoch remained. Recording will be removed.'])
            to_delete{end +1} = EEGtemp.filename;
        end
        
        nEpoch(iRepet) = size(EEGtemp.data,3);
        EEGtempList{iRepet} = EEGtemp;
    end
    
    % extract two extreme cases
    [~, sortInds] = sort(nEpoch);
    maxEp = sortInds(end);
    minEp = sortInds(1);
    selEps = [maxEp, minEp];
    
    numSampEp = 2500; %10s
    %numSampEp = 500; %2s
    sLength = size(EEGtemp.times,2);
    for iRepet = selEps
        figure;
        subLen = 14501;
        selChan = 1;
        plot(EEGselectedBefore{1}.times,EEGselectedBefore{1}.data(selChan,:), 'b');
        
        s210 = contains({EEGtemp0{iRepet}.event.code},'S210');
        latencies210 = [EEGtemp0{iRepet}.event(s210).latency];
        for iEv = 1:sum(s210)
            hold on
            xline(latencies210(iEv)/EEGselectedBefore{1}.srate,'g');
        end
        csm = EEGtemp0{iRepet}.etc.clean_sample_mask;
        hold on
        plot(EEGselectedBefore{1}.times(~csm),EEGselectedBefore{1}.data(1,~csm),'m');
    
        EEGtemp = EEGtempList{iRepet};
        for iEp =1:size(EEGtemp.data,3)
            display(iEp);
            datEp = squeeze(EEGtemp.data(selChan,:, iEp));
            if(iEp ==1)
                searchInds = 1:(size(EEGselectedBefore{1}.data,2)-numSampEp);
            else
                %searchInds = setdiff(searchInds, iMax+1:iMax+1250-1);
                searchInds = iMax+round(numSampEp/2)-2:(size(EEGselectedBefore{1}.data,2)-numSampEp);
            end
            display(size(searchInds))
            corrVal = zeros(1,length(searchInds));
            tic
            for i = searchInds
                ccMat = corrcoef(datEp', squeeze(EEGselectedBefore{1}.data(selChan,i:i+numSampEp-1))');
                corrVal(i) = ccMat(1,2);
                if(corrVal(i)>0.99)
                    fprintf('broke off\n')
                    break
                end
            end
            toc
            [~, iMax] =max(corrVal);
            lat(iEp) = iMax;
            plot(EEGselectedBefore{1}.times(iMax:iMax+sLength-1), datEp+20+2*iEp)
        end
    
%         plot(EEGselectedBefore{1}.data(selChan,1:subLen), 'b');
%         hold on    
%         durTot = 0;
%         for iEv = 1:length(EEGtemp0{iRepet}.event)
%             currLat = EEGtemp0{iRepet}.event(iEv).latency;
%             currDur = EEGtemp0{iRepet}.event(iEv).duration;
%             if(currLat > subLen)
%                 break
%             end
%             if(isnan(currDur))
%                 continue
%             end
%             plot(currLat+durTot:currLat+durTot+currDur-1, EEGselectedBefore{1}.data(selChan,currLat+durTot:currLat+durTot+currDur-1), 'r');
%             durTot = durTot + currDur;
%         end
%     
%         EEGtemp = EEGtempList{iRepet};
%         for iEp =1:size(EEGtemp.data,3)
%             display(iEp)
%             datEp = squeeze(EEGtemp.data(selChan,:, iEp));
%             if(iEp ==1)
%                 searchInds = 1:(size(EEGselectedBefore{1}.data,2)-numSampEp);
%             else
%                 %searchInds = setdiff(searchInds, iMax+1:iMax+1250-1);
%                 searchInds = iMax+round(numSampEp/2)-2:(size(EEGselectedBefore{1}.data,2)-numSampEp);
%             end
%             display(size(searchInds))
%             corrVal = zeros(1,length(searchInds));
%             tic
%             for i = searchInds
%                 ccMat = corrcoef(datEp', squeeze(EEGselectedBefore{1}.data(selChan,i:i+numSampEp-1))');
%                 corrVal(i) = ccMat(1,2);
%                 if(corrVal(i)>0.99)
%                     fprintf('broke off\n')
%                     break
%                 end
%             end
%             toc
%             [~, iMax] =max(corrVal);
%             plot(iMax:iMax+numSampEp-1, datEp+100)
%         end
    end
    
    % Visualize histogram of number of repetitions that resulted in certain
    % number of extracted epochs.
    figure;
    hist(nEpoch)
    xlim([0,max(nEpoch)])
    xlabel('number of epochs')
    ylabel('number of repetitions')
    
    % Save datafile, clear it from memory and store it in the ALLEEG structure
    if exist('clean_channel_mask','var'), EEGtemp.etc.clean_channel_mask = clean_channel_mask; end; clear 'clean_channel_mask';
    EEGtemp = pop_saveset(EEGtemp, 'savemode', 'resave');
    EEGtemp.data = 'in set file';
    ALLEEG = eeg_store(ALLEEG, EEGtemp, iRec);
    STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', params.study, 'filepath', params.preprocessed_data_path);
end
preptime = toc;

% Delete recordings
for iRec = 1:length(to_delete)
    mask = strcmp({ALLEEG.filename},to_delete{iRec});
    ALLEEG(mask) = [];
    STUDY.datasetinfo(mask) = [];
    STUDY.subject(mask) = [];
end
% Save study
STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', [params.study '-clean'], 'filepath', params.preprocessed_data_path);

% PLOTTING PREPROCESSING
% Visualization of detected bad channels. If you set 'FuseChanRej' on, the union of bad channels in all tasks is rejected!
plot_badchannels(params,ALLEEG);
% Visualization of rejected ICs
plot_ICs(params,ALLEEG);
% Visualization of rejected time segments
plot_badtimesegments(params,ALLEEG);

clear EEGtemp ALLEEG;


end

