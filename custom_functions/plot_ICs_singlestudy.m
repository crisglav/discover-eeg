function [IC_plot, ic_kept] = plot_ICs_singlestudy(params,bidsID)
if contains(bidsID,'_')
    x = strsplit(bidsID,'_');
    x = x{1};
else
    x = bidsID;
end
datapath = fullfile(params.PreprocessedDataPath,x,'eeg',[bidsID '_eeg.set']);

% Read header file
hdr = ft_read_header(datapath);
etc = hdr.orig.etc;
nChans = hdr.nChans;

if (isfield(etc,'clean_channel_mask'))
    badchans = ~etc.clean_channel_mask;
else
    badchans = zeros(nChans,1);
end

% Only select the classes that have been marked in params.IClabel (default muscle and eye ICs)
classes = etc.ic_classification.ICLabel.classes;
mask_classes = all(~isnan(params.ICLabel),2);

% IClabel classifications
classifications = etc.ic_classification.ICLabel.orig_classifications;

% Extract the number of bad ICs of each category for all subjects (ICs
% marked bad if they lie between the values of params.IClabel)
badics = sum((classifications(:,2:end) > params.ICLabel(2:end,1)').*(classifications(:,2:end) < params.ICLabel(2:end,2)'));

goodics = nChans - sum(badchans) - sum(badics);
ics = [goodics, badics(mask_classes(2:end)), sum(badchans)];

% Colormap and legend
c = lines(sum(mask_classes)+2); % colormap
c(end,:) = [1 1 1]; % Bad channels to white
l = ['Kept ICs',classes(mask_classes), 'Bad chan']; % Legend

% Deal with the case in which all components that weren't brain were rejected
if all(~isnan(params.ICLabel(1,:)))
    nonbrain_ICs = sum((classifications(:,1) > params.ICLabel(1,1)).*(classifications(:,1) < params.ICLabel(1,2)));
    goodics =  nChans - sum(badchans) - nonbrain_ICs;
    ics = [goodics, nonbrain_ICs, sum(badchans)];
    c = lines(1);
    c(2,:) = [0 0 0]; % Set non-brain components to black
    c(3,:) = [1 1 1]; % Bad channels to white
    l = {'Brain (kept)','Non-brain'};
end

IC_plot =figure('Units','centimeters','Position', [0 0 18 2],'Visible','Off');
h = barh(1,ics,'stacked');
for k = 1:length(c), h(k).FaceColor = c(k,:); end
legend(l,'Location','northeastoutside');
set(gca,'ytick',[],'yticklabel',{''});
box('off')
xlabel('Independent components');
title('IC classification');

ic_kept = goodics/(nChans-sum(badchans));
end

