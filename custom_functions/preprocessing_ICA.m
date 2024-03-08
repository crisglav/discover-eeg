function EEGtemp = preprocessing_ICA(EEGtemp,params)
% Wraper for performing steps 4, 5 and 6

% 4. REMOVE ARTIFACTS WITH ICA
try
    EEGtemp = pop_runica(EEGtemp,'icatype','runica','concatcond','off');
    EEGtemp = pop_iclabel(EEGtemp,'default');
    EEGtemp = pop_icflag(EEGtemp, params.ICLabel); % flag artifactual components using IClabel
    classifications = EEGtemp.etc.ic_classification.ICLabel.classifications; % Keep classifications before component substraction
    EEGtemp = pop_subcomp(EEGtemp,[],0); % Subtract artifactual independent components
    EEGtemp.etc.ic_classification.ICLabel.orig_classifications = classifications;
catch 
    EEGtemp = 4;
    return;
end

% 5. INTERPOLATE MISSING CHANNELS
try
    urchanlocs = EEGtemp.urchanlocs;
    l = length(urchanlocs);
    [~, iref] = setdiff({EEGtemp.chanlocs.labels},{EEGtemp.urchanlocs.labels}); % Handle reference in case it was added back
    if ~isempty(iref)
        urchanlocs(l+1) = EEGtemp.chanlocs(iref);
    end
    EEGtemp = pop_interp(EEGtemp, urchanlocs, 'spherical');
catch 
    EEGtemp = 5;
    return;
end

% 6. REMOVE BAD TIME SEGMENTS
try
    EEGtemp2 = pop_clean_rawdata(EEGtemp,'FlatLineCriterion','off',...
        'ChannelCriterion','off',...
        'LineNoiseCriterion','off',...
        'Highpass','off',...
        'BurstCriterion',params.BurstCriterion,...
        'WindowCriterion',params.WindowCriterion,...
        'BurstRejection','on',...
        'Distance','Euclidian',...
        'WindowCriterionTolerances',params.WindowCriterionTolerances);
    EEGtemp.etc.eventsAfterCRD = EEGtemp.event; % Keep events for visualization later on.
    
    % Either return the data segmented or only the bad segments mask
    if strcmp(params.RejectBadTimeSegments,'off')
        EEGtemp.etc.clean_sample_mask = EEGtemp2.etc.clean_sample_mask; % Return original data and mask of bad time segments
    else
        EEGtemp = EEGtemp2; % Return concatenated data without bad segments and with boundary markers
    end
catch 
    EEGtemp = 6;
    return;
end
end

