function compute_dwpli(params,bidsID,freqBand)
% Load EEG data and check that there are more than 10 epochs
data = load_preprocessed_data(params,bidsID);
assert(size(data.trial,2)>=10, 'Recording with less than 10 epochs');

% Band-pass filter the data in the relevant frequency band
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = params.freq_band.(freqBand);
data = ft_preprocessing(cfg, data);

% Load source reconstruction data
load(fullfile(params.source_folder,[bidsID '_source_' freqBand '.mat']),'source');

% Reconstruct the virtual time series (apply spatial filter to sensor level
% data)
cfg  = [];
cfg.pos = source.pos(source.inside,:);
virtChan_data = ft_virtualchannel(cfg,data,source);
clear data source;

% Frequencies of interest
fois = params.freq_band.(freqBand)(1):params.freq_res_connectivity:params.freq_band.(freqBand)(2);

tic
% Fourier components
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
cfg.pad = 'nextpow2';
cfg.foi = fois;
cfg.tapsmofrq = 1;
virtFreq = ft_freqanalysis(cfg, virtChan_data);
clear virtChan_data;

% Connectivity
cfg = [];
cfg.method = 'wpli_debiased';
source_conn = ft_connectivityanalysis(cfg, virtFreq);
clear virtFreq;

% Average across frequency bins
connMatrix = mean(abs(source_conn.wpli_debiasedspctrm),3);
if all(all(isnan(connMatrix)))
    error('Connectivity matrix only contains NaN');
end
t = toc;
disp([bidsID '_', freqBand, ' computation took ', num2str(t/60), ' minutes'])

save(fullfile(params.connectivity_folder,[bidsID '_dwpli_' freqBand '.mat']),'connMatrix')
end