function plot_badtimesegments(params,EEG)

etc = {EEG.etc};
events = {EEG.event};
nRec = length(EEG);

% Lengths of each recording
lengths = cell2mat(cellfun(@(x) length(x.clean_sample_mask), etc, 'UniformOutput',0)); 

% Keep only boundary events
ev_mask = ~cellfun(@isempty, events);
aux = cellfun(@(x) x(contains({x.type},'boundary')),events(ev_mask),'uni',0); 
events(ev_mask) = aux;

nevents = cell2mat(cellfun(@(x) length(x), events, 'UniformOutput',0));
segs = zeros(nRec,max(nevents)*2+1);

for iSubj=1:nRec
    mask = etc{iSubj}.clean_sample_mask;
    badsegs = reshape(find(diff([false ~mask(:)' false])),2,[]);
    s = [1; badsegs(:); lengths(iSubj)];
    s = s(2:end)-s(1:end-1); 
    segs(iSubj,1:length(s)) = s;   
end

srates = cell2mat({EEG.srate})';
segs = segs./(srates*60); % Divide by sampling rate

f = figure('Position',[1988 548 781 781]);
h = bar(segs,'stacked','EdgeColor','none');
set(h,'FaceColor','Flat');

for iSubj=1:nRec
    for k = 1:find((segs(iSubj,:)~=0),1,'last')
        if mod(k,2) % Deal with the case in which the recording starts with a bad segment
            h(k).CData(iSubj,:) = [0, 0.4470, 0.7410];
        else
            h(k).CData(iSubj,:) = [1,0,0];
        end
    end    
end
yl = ylim;
ylim([0,yl(2)]);

% Recording IDs
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

set(gca,'xticklabel',s_ids,'xticklabelrotation',45)
% xlabel('Recording ID');
ylabel('Length of the recording (minutes)');
legend({'Good','Bad'},'Location','southeast');
title('Segmentation of the data after bad segment removal');
saveas(f,fullfile(params.figures_preprocessing_folder, 'BadSegments.svg'),'svg');

end

