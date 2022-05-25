function [bc_plot] = plot_badchannels(params,bidsID)

% load EEG study info
STUDY = pop_loadstudy('filename', [params.study '_preprocessed.study'], 'filepath', params.preprocessed_data_path);
eeg_idx = find(contains({STUDY.datasetinfo.filename}, bidsID));
EEG = pop_loadset('filepath',STUDY.datasetinfo(eeg_idx).filepath, 'filename',STUDY.datasetinfo(eeg_idx).filename);

etc = {EEG.etc};
nChans = length(EEG(1).urchanlocs);
nRec = 1;

recmask = cellfun(@(x) isfield(x, 'clean_channel_mask'), etc); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_channel_mask, etc(recmask),'uni',0);
badrecs = double(~cat(2,tmp{:}));
badchans = zeros(nChans,nRec);
badchans(:,recmask) = badrecs;

f = figure('Position',[1988 548 781 100], 'visible', 'off');

% Channel labels
chanlabels = {EEG(1).urchanlocs.labels};

% colormap
cmap = [1 1 1; lines(1)];
if all(all(badchans == 0)), cmap = [1 1 1]; end

% Plot the heatmap
bc_plot = heatmap(chanlabels, '-', badchans.','CellLabelColor','none','Colormap',cmap,'ColorbarVisible','off');
bc_plot.Title = 'Bad channels';
% h.XLabel = 'Recording ID';
bc_plot.XLabel = 'Channels';

end

