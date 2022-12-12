function compute_aec(params,bidsID,freqBand)
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

% Compute AEC
s = tic;
connMatrix = aecConnectivity(virtChan_data);
connMatrix = mean(connMatrix,3);
t = toc(s);
disp([bidsID '_', freqBand, ' computation took ', num2str(t/60), ' minutes'])

save(fullfile(params.connectivity_folder,[bidsID '_aec_' freqBand '.mat']),'connMatrix')

end