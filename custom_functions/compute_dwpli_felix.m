function compute_dwpli_felix(params,bidsID,freqBand)
% Instead of computing one connectivity matrix per frequency bin, obtain
% the fourier time series with multitapering with a high frequency
% resolution and obtain the connectivity matrix at that freq.

% Load EEG data
data = load_preprocessed_data(params,bidsID);

% Load source reconstruction data
load(fullfile(params.source_folder,[bidsID '_source_' freqBand '.mat']),'source');

% Reconstruct the virtual time series (apply spatial filter to sensor level
% data)
nElec =  numel(data.label);
nVoxel = length([source.avg.filter{:}])/nElec;

filter_array = reshape([source.avg.filter{:}], nElec, nVoxel)';
virtChan_timeSeries = cellfun(@(x) filter_array*x, data.trial, 'UniformOutput', false);

% convert the reconstructed time series into correct format
virtChan_data = data;
virtChan_data.label = cellstr(num2str(source.pos(source.inside,:)));
virtChan_data.trial = virtChan_timeSeries;
clear virtChan_timeSeries sourceData filter_array;

% prealocate connectivity matrix
connMatrix = zeros(nVoxel, nVoxel);

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

% Connectivity
cfg = [];
cfg.method = 'wpli_debiased';
source_conn = ft_connectivityanalysis(cfg, virtFreq);
connMatrix(:,:) = source_conn.(['wpli_debiased','spctrm']);
b = toc;

[nRows, nCols] = size(connMatrix);
connMatrix(1:nRows+1:nRows*nCols) = nan;
disp([bidsID '_', freqBand, ' computation took ', num2str(b/60), ' minutes'])
save(fullfile(params.connectivity_folder,[bidsID '_dwpliFelix_' freqBand '.mat']),'connMatrix')

end