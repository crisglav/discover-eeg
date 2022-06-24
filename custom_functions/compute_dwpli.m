function connMatrix = compute_dwpli(params,bidsID,freqBand)
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

% Frequencies of interest
fois = params.freq_band.(freqBand)(1):params.freq_res_connectivity:params.freq_band.(freqBand)(2);

% prealocate connectivity matrix
connMatrix = zeros(nVoxel, nVoxel, length(fois));

tic
for idx=1:length(fois)
    % Fourier components
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.taper = 'dpss';
    cfg.output = 'fourier';
    cfg.keeptrials = 'yes';
    cfg.pad = 'nextpow2';
    cfg.foi = fois(idx);
    cfg.tapsmofrq = 1;
    virtFreq = ft_freqanalysis(cfg, virtChan_data);

    % Connectivity
    cfg = [];
    cfg.method = 'wpli_debiased';
    source_conn = ft_connectivityanalysis(cfg, virtFreq);
    connMatrix(:,:,idx) = source_conn.(['wpli_debiased','spctrm']);
end
t = toc;

connMatrix = mean(abs(connMatrix), 3);
[nRows, nCols] = size(connMatrix);
connMatrix(1:nRows+1:nRows*nCols) = nan;
disp([bidsID '_', freqBand, ' computation took ', num2str(t/60), ' minutes'])
clear connMatrix;

save(fullfile(params.connectivity_folder,[bidsID '_dwpli_' freqBand '.mat']),'connMatrix')
end