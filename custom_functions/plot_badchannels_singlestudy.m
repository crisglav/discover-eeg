function [bc_plot] = plot_badchannels_singlestudy(params,bidsID)
x = strsplit(bidsID,'_');
x = x(1:end-1);
datapath = fullfile(params.preprocessed_data_path,x{:},'eeg',[bidsID '_eeg.set']);

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
cmap = [1 1 1; lines(1)];
if all(all(badchans == 0)), cmap = [1 1 1]; end

% Plot the heatmap
figure('Position',[1988 548 1500 200],'visible','off');
bc_plot = heatmap(chanlabels,'-',double(badchans)','CellLabelColor','none','Colormap',cmap,'ColorbarVisible','off');
bc_plot.Title = 'Bad channels';
bc_plot.XLabel = 'Channels';
end

