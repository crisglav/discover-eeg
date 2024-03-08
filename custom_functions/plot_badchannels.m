function plot_badchannels(params,EEG)

% Recording ids
s_ids = cellfun(@(x) regexp(x,'.*(?=_eeg.set)','match','lineanchors'),{EEG.filename});

% Check that the number of channels is the same for all the recordings
etc = {EEG.etc};
if range([EEG.nbchan]) ~= 0
    error('Different number of total channels in at least one recording')
end
nChans = length(EEG(1).urchanlocs);
chanlabels = {EEG(1).urchanlocs.labels};  % Channel labels

% Plot in batches of x recordings
x = 30;
nRec = length(EEG);
nBatches = ceil(nRec/x);

for iBatch=1:nBatches
    
    % Subset of recordings for this batch
    if nRec > x*iBatch
        bmask = (x*(iBatch-1)+1):x*iBatch;
    else
        bmask = (x*(iBatch-1)+1):nRec;
    end
    
    % Recording ids for this batch
    % For a short version of the recordings IDs: delete fields that are all the
    % same (e.g. if all sessions are ses-1)
    splitted_ids = cellfun(@(x) strsplit(x,'_'),s_ids(bmask),'UniformOutput',false);
    splitted_ids = vertcat(splitted_ids{:});
    mask = ones(1,size(splitted_ids,2));
    for i=1:size(splitted_ids,2)
        if (numel(unique(splitted_ids(:,i)))==1 && size(splitted_ids,1)>1)
            mask(i)=0;
        end
    end
    ids = splitted_ids(:,find(mask));
    ids = join(ids,'_',2);
    ids = insertBefore(ids,'_','\'); % Escape the underscores
    
    
    % Bad channels for this batch
    % Extract the clean_channel_mask for all recordings
    betc = etc(bmask);
    recmask = cellfun(@(x) isfield(x, 'clean_channel_mask'), betc); % Deal with the case where no bad channels were detected
    tmp = cellfun(@(x) x.clean_channel_mask, betc(recmask),'uni',0);
    tmp = cellfun(@(x) x(:), tmp, 'UniformOutput', false); % force to be a column vector
    badrecs = double(~cat(2,tmp{:}));

    badchans = zeros(nChans,length(bmask));
    badchans(:,recmask) = badrecs;

    % colormap
    cmap = [1 1 1; 1 0 0];
    if all(all(badchans == 0)), cmap = [1 1 1]; end

    % Plot the heatmap
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    h = heatmap(ids,chanlabels,badchans,'CellLabelColor','none','Colormap',cmap,'ColorbarVisible','off');
    h.Title = 'Bad channels';
    h.YLabel = 'Channels';
    
    saveas(f,fullfile(params.FiguresPreprocessingPath, ['badchans_' num2str(iBatch) '.svg']),'svg');
    save(fullfile(params.FiguresPreprocessingPath, ['badchans_' num2str(iBatch) '.mat']),'badchans','ids','chanlabels');
end

end

