function power = compute_power(params,bidsID)
% Load EEG data
data = load_preprocessed_data(params,bidsID);

% Average power spectrum across epochs in the range 1-100 Hz
cfg = [];
cfg.foi = 1:params.freq_res:100;
cfg.method = 'mtmfft';
cfg.taper = params.taper;
cfg.tapsmofrq = params.tapsmofrq;
cfg.output = 'pow';
cfg.keeptrials ='no';
power = ft_freqanalysis(cfg, data);
end