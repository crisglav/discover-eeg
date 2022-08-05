function compute_dwpli_felix(params,bidsID,freqBand)
% Instead of computing one connectivity matrix per frequency bin, obtain
% the fourier time series with multitapering with at the middle frequency
% point of the band with a smoothing equal to half of the band range.
% Comment: it makes connectiviy values higher than the other method and
% some values present a bias towards negative values, therefore there are
% negative values in the connectivity matrix.
% Cristina Gil and Felix Bott, 05.08.2022, TUM

% Load EEG data
data = load_preprocessed_data(params,bidsID);

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