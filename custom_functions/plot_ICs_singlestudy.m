function [IC_plot, ic_kept, ic_kept_pc] = plot_ICs_singlestudy(params,bidsID)

% load EEG study info
STUDY = pop_loadstudy('filename', [params.study '_preprocessed.study'], 'filepath', params.preprocessed_data_path);
eeg_idx = find(contains({STUDY.datasetinfo.filename}, bidsID));
EEG = pop_loadset('filepath',STUDY.datasetinfo(eeg_idx).filepath, 'filename',STUDY.datasetinfo(eeg_idx).filename);

etc = {EEG.etc};
nChans = length(EEG.urchanlocs);
nRec = 1;

% Bad channels
recmask = cellfun(@(x) isfield(x, 'clean_channel_mask'), etc); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_channel_mask, etc(recmask),'uni',0);
badrecs = double(~cat(2,tmp{:}));
badchans = zeros(nChans,nRec);
badchans(:,recmask) = badrecs;

% Only select the classes that have been marked in params.IClabel (default muscle and eye ICs)
classes = etc{1}.ic_classification.ICLabel.classes;
mask_classes = all(~isnan(params.IClabel),2);

% IClabel classifications
classifications = cellfun(@(x) x.ic_classification.ICLabel.orig_classifications, etc,'uni',0);
[~, classifications_array] = max(classifications{1}, [], 2);
% Count rejected ICs according to the pipeline parameters
extract_bad_ICs = @(x) sum((x(:,2:end) > params.IClabel(2:end,1)').*(x(:,2:end) < params.IClabel(2:end,2)'));
badics = cellfun(extract_bad_ICs, classifications,'uni',0);
badics = reshape(cell2mat(badics),length(classes)-1,nRec);
% X data
totalbadics = sum(badics,1);
goodics = nChans - sum(badchans) - totalbadics;
ics = [goodics;badics(mask_classes(2:end),:);sum(badchans)]';
ic_kept = goodics;
ic_kept_pc = goodics / (sum(ics) - sum(badchans));
% % Colormap and legend
c = lines(sum(mask_classes)+2); % colormap
c(end,:) = [1 1 1]; % Bad channels to white
l = ['Kept ICs',classes(mask_classes), 'Bad channels']; % Legend
if all(~isnan(params.IClabel(1,:)))
    extract_nonbrain_ICs = @(x) sum((x(:,1) > params.IClabel(1,1)).*(x(:,1) < params.IClabel(1,2)));
    extract_brain_ICs = @(x) sum(~((x(:,1) > params.IClabel(1,1)).*(x(:,1) < params.IClabel(1,2))));
    badics =  cellfun(extract_nonbrain_ICs, classifications);
    goodics =  cellfun(extract_brain_ICs, classifications);
    ics = [goodics;badics;sum(badchans)]';
    c(end-1,:) = [0 0 0]; % Set non-brain components to black
    l = {'Brain (kept)','Non-brain'};
end


% Create stacked bar plot
IC_plot = figure('Position',[1988 548 781 781],'visible', 'off');
p = pie(ics);
IC_plot.Colormap = c;
T = p(strcmpi(get(p, 'Type'), 'text'));
P = cell2mat(get(T, 'Position'));
set(T, {'Position'}, num2cell(P*0.5,2));
text(P(:,1), P(:,2), strcat(l, ': ', string(ics)), 'BackgroundColor', 'white');
ylabel('Independent components');
title('IC classification');

end

