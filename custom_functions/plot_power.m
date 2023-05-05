function [power_fig, topoplot_fig] = plot_power(params,bidsID)

% Load pre-computed power and peak frequency
load(fullfile(params.PowerPath,[bidsID '_power.mat']),'power');
load(fullfile(params.PowerPath,[bidsID '_peakfrequency.mat']),'peakfrequency')

% Average power across channels
avgpow = mean(power.powspctrm,1);

freqNames = fields(params.FreqBand)';    
% Plot average power spectrum across channels
power_fig = figure('Units','centimeters','Position', [0 0 12 7],'Visible','off');
ax1 = axes(power_fig);
plot(ax1, power.freq,avgpow,'k');
hold on
% Color the different frequency bands
c = lines(length(freqNames));
for iFreq =1:length(freqNames)
    freqRange = find(power.freq >= params.FreqBand.(freqNames{iFreq})(1) & power.freq <= params.FreqBand.(freqNames{iFreq})(2));
    a(iFreq) = area(power.freq(freqRange),avgpow(freqRange),'FaceColor',c(iFreq,:),'DisplayName',freqNames{iFreq}); % Color the area
    hold on;
end
legend(a)

% Plot APF over the average power spectrum
ax2 = copyobj(ax1,gcf);
if ~isempty(peakfrequency.localmax)
    p(1) = plot(ax2,peakfrequency.localmax,interp1(power.freq,avgpow,peakfrequency.localmax),'v','MarkerSize',6,'MarkerEdgeColor','#006633','MarkerFaceColor','#006633','DisplayName','APF - max');
    hold on
    p(2) =  plot(ax2,peakfrequency.cog,interp1(power.freq,avgpow,peakfrequency.cog),'v','MarkerSize',6,'MarkerEdgeColor','#95C11F','MarkerFaceColor','#95C11F','DisplayName','APF - c.o.g.');
else
    p = plot(ax2,peakfrequency.cog,interp1(power.freq,avgpow,peakfrequency.cog),'v','MarkerSize',6,'MarkerEdgeColor','#95C11F','MarkerFaceColor','#95C11F','DisplayName','APF - c.o.g.');
end
legend(p,'Location','east','Color','none')

title('Power spectrum (electrode avg.)');
ylabel('Power (uV^2/Hz)');
xlabel('Frequency (Hz)');
box off
set(ax2,'Color','none','box','off','Position',ax1.Position);

% Topoplots of different frequency bands
topoplot_fig = figure('Units','centimeters','Position', [0 0 11 10], 'Visible','Off');
tcl = tiledlayout(2,2);
tcl.TileSpacing = 'compact';
tcl.Padding = 'compact';
for iFreq=1:length(freqNames)
    cfg = [];
    cfg.xlim = params.FreqBand.(freqNames{iFreq});
%     cfg.marker = 'labels'; % Uncomment if you want the electrode names
    cfg.comment = 'no';
    ft_topoplotER(cfg,power)
    title(freqNames{iFreq});
    cb = colorbar;
    ax = gca;
    ax.Parent = tcl;
    ax.Layout.Tile = iFreq;
    cb.Parent = tcl;
    close;
end

end