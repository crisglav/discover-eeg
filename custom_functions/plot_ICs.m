function plot_ICs(params,EEG)
% Recording ids
s_ids = cellfun(@(x) regexp(x,'.*(?=_eeg.set)','match','lineanchors'),{EEG.filename});

% Check that the number of channels is the same for all the recordings
etc = {EEG.etc};
if range([EEG.nbchan]) ~= 0
    error('Different number of total channels in at least one recording')
end
nChans = length(EEG(1).urchanlocs);

% Only select the IC classes that have been marked in params.IClabel (default muscle and eye ICs)
classes = etc{1}.ic_classification.ICLabel.classes;
mask_classes = all(~isnan(params.ICLabel),2);
extract_bad_ICs = @(x) sum((x(:,2:end) > params.ICLabel(2:end,1)').*(x(:,2:end) < params.ICLabel(2:end,2)'));

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
    
    % IClabel classifications
    classifications = cellfun(@(x) x.ic_classification.ICLabel.orig_classifications, betc,'uni',0);
    
    % Extract the number of bad ICs of each category for all subjects (ICs
    % marked bad if they lie between the values of params.IClabel)
    badics = cellfun(extract_bad_ICs, classifications,'uni',0);
    badics = reshape(cell2mat(badics),length(classes)-1,length(bmask));
    
    % create a structure with the good ICs, badIcs and bad channels
    totalbadics = sum(badics,1);
    goodics = nChans - sum(badchans) - totalbadics;
    ics = [goodics;badics(mask_classes(2:end),:);sum(badchans)]';
    
    % Colormap and legend
    c = lines(sum(mask_classes)+2); % colormap
    c(end,:) = [1 1 1]; % Bad channels to white
    l = ['Kept ICs',classes(mask_classes)]; % Legend
    
    % Deal with the case in which all components that weren't brain were rejected
    if all(~isnan(params.ICLabel(1,:)))
        extract_nonbrain_ICs = @(x) sum((x(:,1) > params.ICLabel(1,1)).*(x(:,1) < params.ICLabel(1,2)));
        extract_brain_ICs = @(x) sum(~((x(:,1) > params.ICLabel(1,1)).*(x(:,1) < params.ICLabel(1,2))));
        badics =  cellfun(extract_nonbrain_ICs, classifications);
        goodics =  cellfun(extract_brain_ICs, classifications);
        ics = [goodics;badics;sum(badchans)]';
        c(end-1,:) = [0 0 0]; % Set non-brain components to black
        l = {'Brain (kept)','Non-brain'};
    end
    
    % Create stacked bar plot
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    h = bar(1:size(ics,1), ics,'stacked','EdgeColor','none');
    for k = 1:length(c), h(k).FaceColor = c(k,:); end
    xticks(1:length(bmask));
    box off;
    set(gca,'xticklabel',ids(1:length(bmask)),'xticklabelrotation',45);
    ylabel('Independent components');
    legend(l,'Location','southeast');
    title('IC classification');
    
    saveas(f,fullfile(params.FiguresPreprocessingPath, ['ICA_' num2str(iBatch) '.svg']),'svg');
    save(fullfile(params.FiguresPreprocessingPath, ['ICA_' num2str(iBatch) '.mat']),'ics','ids','l');
end

end

