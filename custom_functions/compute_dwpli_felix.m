function compute_dwpli_felix(params,bidsID,freqBand)
% Instead of computing one connectivity matrix per frequency bin, obtain
% the fourier time series with multitapering with a high frequency
% resolution and obtain the connectivity matrix at that freq.

% Load EEG data
data = load_preprocessed_data(params,bidsID);

% Band-pass filter the data in the relevant frequency band
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = params.freq_band.(freqBand);
data = ft_preprocessing(cfg, data);

% Load source reconstruction data
load(fullfile(params.source_folder,[bidsID '_source_' freqBand '.mat']),'source');
nVoxel = sum(source.inside);

% Reconstruct the virtual time series (apply spatial filter to sensor level
% data)
cfg  = [];
cfg.pos = source.pos(source.inside,:);
virtChan_data = ft_virtualchannel(cfg,data,source);
clear data source;

% Frequency that is in the middle of the frequency band of interest
halfrange = (params.freq_band.(freqBand)(2)-params.freq_band.(freqBand)(1))/2;
midfreq = params.freq_band.(freqBand)(1) + halfrange;

tic
% Fourier components
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'dpss'; 
cfg.tapsmofrq = halfrange; % set the spectral smoothing to half of the frequency band range
cfg.foi = midfreq;
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
cfg.pad = 'nextpow2';
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

t = toc;
disp([bidsID '_', freqBand, ' computation took ', num2str(t/60), ' minutes'])

save(fullfile(params.connectivity_folder,[bidsID '_dwpli-Felix_' freqBand '.mat']),'connMatrix')
end