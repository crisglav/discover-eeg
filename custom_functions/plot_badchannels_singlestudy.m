function [bc_plot] = plot_badchannels_singlestudy(params,bidsID)
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

% Channel labels
chanlabels = {hdr.orig.urchanlocs.labels};

% colormap
cmap = [1 1 1; 1 0 0];
if all(all(badchans == 0)), cmap = [1 1 1]; end

% Plot the heatmap
figure('Units','centimeters','Position', [0 0 18 2],'Visible','Off');
bc_plot = heatmap(chanlabels,'-',double(badchans)','CellLabelColor','none','Colormap',cmap,'ColorbarVisible','off');
bc_plot.Title = 'Bad channels';
bc_plot.XLabel = 'Channel names';
end

