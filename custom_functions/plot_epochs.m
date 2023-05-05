function plot_epochs(params,EEG)
nepochs = cellfun(@length,{EEG.epoch});
nRec = length(EEG);

f = figure('units','normalized','outerposition',[0 0 1 1]);
bar(nepochs);

% Recording IDs
s_ids = cellfun(@(x) regexp(x,'.*(?=_eeg.set)','match','lineanchors'),{EEG.filename});
% For a short version of the recordings IDs: delete fields that are all the
% same (e.g. if all sessions are ses-1)
splitted_ids = cellfun(@(x) strsplit(x,'_'),s_ids,'UniformOutput',false);
splitted_ids = vertcat(splitted_ids{:});
mask = ones(1,size(splitted_ids,2));
for i=1:size(splitted_ids,2)
    if (numel(unique(splitted_ids(:,i)))==1 && size(splitted_ids,1)>1)
        mask(i)=0;
    end
end
s_ids = splitted_ids(:,find(mask));
s_ids = join(s_ids,'_',2);
s_ids = insertBefore(s_ids,'_','\'); % Escape the underscores

xticks(1:2:nRec);
box off;
set(gca,'xticklabel',s_ids(1:2:nRec),'xticklabelrotation',45)
xlabel('Recording ID');
ylabel('Number of clean epochs');
title('Number of epochs after preprocessing');
saveas(f,fullfile(params.FiguresPreprocessingPath, 'Epochs.svg'),'svg');
% savefig(f,fullfile(params.FiguresPreprocessingPath, 'Epochs.fig'));
save(fullfile(params.FiguresPreprocessingPath, 'Epochs.mat'),'nepochs');

end