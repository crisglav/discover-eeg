function plot_badchannels(params,EEG)

etc = {EEG.etc};
nChans = length(EEG(1).urchanlocs);
nRec = length(EEG);

recmask = cellfun(@(x) isfield(x, 'clean_channel_mask'), etc); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_channel_mask, etc(recmask),'uni',0);
badrecs = double(~cat(2,tmp{:})); % error here when you have datasets with different number of electrodes
badchans = zeros(nChans,nRec);
badchans(:,recmask) = badrecs;

f = figure('Position',[1988 548 781 781], 'visible', 'off');

% Recording ids
s_ids = cellfun(@(x) regexp(x,'.*(?=_eeg.set)','match','lineanchors'),{EEG.filename});
% For a short version of the recordings IDs: delete fields that are all the
% same (e.g. if all sessions are ses-1)
splitted_ids = cellfun(@(x) strsplit(x,'_'),s_ids,'UniformOutput',false);
splitted_ids = vertcat(splitted_ids{:});
mask = ones(1,size(splitted_ids,2));
for i=1:size(splitted_ids,2)
    if numel(unique(splitted_ids(:,i)))==1
        mask(i)=0;
    end
end
s_ids = splitted_ids(:,find(mask));
s_ids = join(s_ids,'_',2);
s_ids = insertBefore(s_ids,'_','\'); % Escape the underscores

% Channel labels
chanlabels = {EEG(1).urchanlocs.labels};

% colormap
cmap = [1 1 1; lines(1)];
if all(all(badchans == 0)), cmap = [1 1 1]; end

% Plot the heatmap
h = heatmap(s_ids,chanlabels,badchans,'CellLabelColor','none','Colormap',cmap,'ColorbarVisible','off');
h.Title = 'Bad channels';
% h.XLabel = 'Recording ID';
h.YLabel = 'Channels';
saveas(f,fullfile(params.figures_preprocessing_folder, 'badchans.svg'),'svg');

end

