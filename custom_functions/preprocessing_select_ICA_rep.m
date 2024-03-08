function EEGtemp = preprocessing_select_ICA_rep(EEGtemp_clean)

% Select the repetition closer to the 'average' bad time segments mask
recmask = cellfun(@(x) isfield(x, 'etc'), EEGtemp_clean);
etc = cellfun(@(x) x.etc, EEGtemp_clean(recmask), 'UniformOutput',0);
csm = cell2mat(cellfun(@(x) x.clean_sample_mask, etc, 'UniformOutput',0)');
csmAverage = mean(csm,1);
distToAvg = sum(abs(csm - csmAverage), 2);
[~,selRun] = min(distToAvg);
ix = find(recmask);
EEGtemp = EEGtemp_clean{ix(selRun)};
           
end

