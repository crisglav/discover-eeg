function compute_peakfrequency(params,bidsID)
% Computes peak frequency with two methods: local maximum and center of
% gravity
% Load pre-computed power spectrum
load(fullfile(params.PowerPath,[bidsID '_power.mat']),'power')

% Average power across channels
avgpow = mean(power.powspctrm,1);

% Frequency range (search limits for the peak = alpha band)
freqRange = find(power.freq >= params.FreqBand.alpha(1) & power.freq <= params.FreqBand.alpha(2));

% Peak frequency computed on the power spectrum averaged across channels
% Approach 1. Find highest local maximum in the power spectrum averaged across epochs
[~, pf_localmax] = findpeaks(avgpow(freqRange),power.freq(freqRange),'SortStr','descend','NPeaks',1);
peakfrequency.localmax = pf_localmax;

% Approach 2. Compute center of gravity on the power spectrum averaged across epochs
pf_cog = sum(avgpow(freqRange).*power.freq(freqRange))/sum(avgpow(freqRange));
peakfrequency.cog = pf_cog;

save(fullfile(params.PowerPath,[bidsID '_peakfrequency.mat']),'peakfrequency')

% % Plotting
% pf_fig = plot_peakfrequency(params,bidsID);
% saveas(pf_fig,fullfile(params.PowerPath,[bidsID '_peakfrequency.svg']));
% close(pf_fig)

% Plotting power
[power_fig, topoplot_fig] = plot_power(params,bidsID);
saveas(power_fig,fullfile(params.PowerPath,[bidsID '_power.svg']));
saveas(topoplot_fig,fullfile(params.PowerPath,[bidsID '_power_topoplots.svg']));
end