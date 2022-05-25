function pf = compute_peak_frequency(params,bidsID)
% Computes peak frequency with two methods: local maximum and center of
% gravity
% Load pre-computed power spectrum
load(fullfile(params.power_folder,[bidsID '_power.mat']),'power')

% Average power across channels
avgpow = mean(power.powspctrm,1);

% Frequency range (search limits for the peak = alpha band)
freqRange = find(power.freq >= params.freq_band.alpha(1) & power.freq <= params.freq_band.alpha(2));

% Peak frequency computed on the power spectrum averaged across channels
% Approach 1. Find highest local maximum in the power spectrum averaged across epochs
[~, pf_localmax] = findpeaks(avgpow(freqRange),power.freq(freqRange),'SortStr','descend','NPeaks',1);
pf.localmax = pf_localmax;

% Approach 2. Compute center of gravity on the power spectrum averaged across epochs
pf_cog = sum(avgpow(freqRange).*power.freq(freqRange))/sum(avgpow(freqRange));
pf.cog = pf_cog;
end