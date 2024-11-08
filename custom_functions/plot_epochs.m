function plot_epochs(params,EEG)

% Recording IDs
s_ids = cellfun(@(x) regexp(x,'.*(?=_eeg.set)','match','lineanchors'),{EEG.filename});

nEpochs = cellfun(@length,{EEG.epoch});

% Plot in batches of x recordings
x = 30;
nRec = length(EEG);
nBatches = ceil(nRec/x);
clear EEG;
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
    
    % Bar plot
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    bar(nEpochs(bmask));
    xticks(1:length(bmask));
    box off;
    set(gca,'xticklabel',ids(1:length(bmask)),'xticklabelrotation',45)
    ylabel('Number of clean epochs');
    title('Number of epochs after preprocessing');
    saveas(f,fullfile(params.FiguresPreprocessingPath, ['Epochs_ ' num2str(iBatch) '.svg']),'svg');
    nEpochs_batch = nEpochs(bmask);
    save(fullfile(params.FiguresPreprocessingPath, ['Epochs_' num2str(iBatch) '.mat']),'nEpochs_batch');
end
end