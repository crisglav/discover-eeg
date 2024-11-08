function plot_badtimesegments(params,EEG)
% Recording ids
s_ids = cellfun(@(x) regexp(x,'.*(?=_eeg.set)','match','lineanchors'),{EEG.filename});

etc = {EEG.etc};
srates = {EEG.srate};

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
    orig_ids = s_ids(bmask);
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
    
    
    % Lengths of each recording
    betc = etc(bmask);        
    lengths = cell2mat(cellfun(@(x) length(x.clean_sample_mask), betc, 'UniformOutput',0)); 
    
    badsegs_batch = cell(1,length(bmask));
    for iRec = 1:length(bmask)
        mask = betc{iRec}.clean_sample_mask;
        % Find transitions 0 to 1 and 1 to 0
        boundaries = find(diff([false ~mask(:)' false]));
        % Add first sample and last sample
        boundaries = [1, boundaries, lengths(iRec)];
        % Substract subsequent boundaries to get the duration of bad time
        % periods
        badsegs_iRec = boundaries(2:end) - boundaries(1:end-1);
        % Store bad timeperiods of this recording in a cell array
        badsegs_batch{iRec} = badsegs_iRec;
        
        
    end
    
    % Find the maximum number of bad segments per recording in this batch
    % and create a matrix with so many rows
    nbadsegs = cellfun(@(x) length(x), badsegs_batch);
    segs = zeros(length(bmask),max(nbadsegs));
    for iRec = 1:length(bmask)
        segs(iRec,1:nbadsegs(iRec)) = badsegs_batch{iRec};
    end
    
    % Divide by sampling rate
    bsrates = cell2mat(srates(bmask))';
    segs = segs./(bsrates*60); % Divide by sampling rate    
    
    % Bar plots
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    h = bar(1:size(segs,1), segs,'stacked','EdgeColor','none');
    set(h,'FaceColor','Flat');
    
    for iRec=1:length(bmask)
        for k = 1:find((segs(iRec,:)~=0),1,'last')
            if mod(k,2) % Deal with the case in which the recording starts with a bad segment
                h(k).CData(iRec,:) = [0, 0.4470, 0.7410];
            else
                h(k).CData(iRec,:) = [1,0,0];
            end
        end
    end
    yl = ylim;
    ylim([0,yl(2)]);
    xticks(1:length(bmask));
    box off;
    set(gca,'xticklabel',ids(1:length(bmask)),'xticklabelrotation',45)
    ylabel('Length of the recording (minutes)');
    legend({'Good','Bad'},'Location','southeast');
    if strcmp(params.RejectBadTimeSegments,"on")
        title('Rejected bad segments');
    else
        title('Detected bad segments (not rejected)');
    end

    saveas(f,fullfile(params.FiguresPreprocessingPath, ['BadSegments_' num2str(iBatch) '.svg']),'svg');
    save(fullfile(params.FiguresPreprocessingPath, ['BadSegments_' num2str(iBatch) '.mat']),'segs','orig_ids');
    close(f);
end


end

