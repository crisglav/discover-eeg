function pf_fig = plot_peakfrequency(params,bidsID)
% Load pre-computed power spectrum and peak frequency
load(fullfile(params.power_folder,[bidsID '_power.mat']),'power')
load(fullfile(params.pf_folder,[bidsID '_peakfrequency.mat']),'pf')

% Average power across channels
avgpow = mean(power.powspctrm,1);

% Frequency range (search limits for the peak = alpha band)
freqRange = find(power.freq >= params.freq_band.alpha(1) & power.freq <= params.freq_band.alpha(2));

% Plotting peak frequency

pf_fig = figure('Position',[1988 548 781 781]);
findpeaks(avgpow(freqRange),power.freq(freqRange),'SortStr','descend','NPeaks',1);
hold on
plot(pf.cog,interp1(power.freq(freqRange),avgpow(freqRange),pf.cog),'v','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319');
title(['Peak frequency - ' bidsID],'Interpreter','None');
ylabel('Power (uV/Hz)')
xlabel('Frequency (Hz)')
hold on

% Hard code the legend
h(1) = plot(nan,nan,'v','MarkerEdgeColor','#0072BD','MarkerFaceColor','#0072BD');
h(2) = plot(nan,nan,'v','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319');
legend(h,{'Maximum peak','Center of gravity'},'Location','southeast');
hold off


end