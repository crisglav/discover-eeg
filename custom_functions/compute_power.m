function compute_power(params,bidsID)
% Load EEG data
data = load_preprocessed_data(params,bidsID);

% Average power spectrum across epochs in the range 1-100 Hz
cfg = [];
cfg.foilim = [1 100];
cfg.method = 'mtmfft';
cfg.taper = params.Taper;
cfg.tapsmofrq = params.Tapsmofrq;
if ~isempty(params.Pad)
    cfg.pad = params.Padding;
    cfg.padtype = 'zero';
end
cfg.output = 'pow';
cfg.keeptrials ='no';
power = ft_freqanalysis(cfg, data);

save(fullfile(params.PowerPath,[bidsID '_power.mat']),'power')
end