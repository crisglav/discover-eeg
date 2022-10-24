function pf_fig = plot_peakfrequency(params,bidsID)
% Load pre-computed power spectrum and peak frequency
load(fullfile(params.power_folder,[bidsID '_power.mat']),'power')
load(fullfile(params.power_folder,[bidsID '_peakfrequency.mat']),'peakfrequency')

% Average power across channels
avgpow = mean(power.powspctrm,1);

% Frequency range (search limits for the peak = alpha band)
freqRange = find(power.freq >= params.freq_band.alpha(1) & power.freq <= params.freq_band.alpha(2));

% Plotting peak frequency
pf_fig = figure('Position',[1988 548 781 781], 'visible', 'off');
plot(power.freq(freqRange),avgpow(freqRange),'k');
hold on
h(1) = plot(peakfrequency.localmax,interp1(power.freq(freqRange),avgpow(freqRange),peakfrequency.localmax),'v','MarkerSize',8,'MarkerEdgeColor','#006633','MarkerFaceColor','#006633');
hold on
h(2) = plot(peakfrequency.cog,interp1(power.freq(freqRange),avgpow(freqRange),peakfrequency.cog),'v','MarkerSize',8,'MarkerEdgeColor','#95C11F','MarkerFaceColor','#95C11F');
ylabel('Power (uV/Hz)')
xlabel('Frequency (Hz)')
legend(h,{'Maximum peak','Center of gravity'},'Location','southeast');
box off
hold off

title(['Peak frequency of ' bidsID],'Interpreter','None');

end