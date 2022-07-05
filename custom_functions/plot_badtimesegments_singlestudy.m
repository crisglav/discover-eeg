function [bs_plot, bs_secs, bs_pc] = plot_badtimesegments_singlestudy(params,bidsID)
x = strsplit(bidsID,'_');
x = x(1:end-1);
datapath = fullfile(params.preprocessed_data_path,x{:},'eeg',[bidsID '_eeg.set']);

hdr = ft_read_header(datapath);
etc = hdr.orig.etc;
events = hdr.orig.event;

if (isfield(etc,'clean_sample_mask'))
    lengths = length(etc.clean_sample_mask);
else
    lengths = zeros(nChans,1);
end

% load EEG study info
STUDY = pop_loadstudy('filename', [params.study '_preprocessed.study'], 'filepath', params.preprocessed_data_path);
eeg_idx = find(contains({STUDY.datasetinfo.filename}, bidsID));
EEG = pop_loadset('filepath',STUDY.datasetinfo(eeg_idx).filepath, 'filename',STUDY.datasetinfo(eeg_idx).filename);

etc = {EEG.etc};
events = {EEG.event};
nRec = 1;

% Lengths of each recording
lengths = cell2mat(cellfun(@(x) length(x.clean_sample_mask), etc, 'UniformOutput',0)); 

% Keep only boundary events
ev_mask = ~cellfun(@isempty, events);
aux = cellfun(@(x) x(contains({x.type},'boundary')),events(ev_mask),'uni',0); 
events(ev_mask) = aux;

nevents = cell2mat(cellfun(@(x) length(x), events, 'UniformOutput',0));
segs = zeros(2,max(nevents)*2+1);

mask = etc{1}.clean_sample_mask;
badsegs = reshape(find(diff([false ~mask(:)' false])),2,[]);
s = [1; badsegs(:); lengths(1)];
s = s(2:end)-s(1:end-1); 
segs(1,1:length(s)) = s;

srates = cell2mat({EEG.srate})';
segs = segs./(srates*60); % Divide by sampling rate

bs_plot = figure('Position',[1988 548 1500 300], 'visible', 'off');
b = barh(1, segs(1,:),'stacked','EdgeColor','none');
set(b,'FaceColor','Flat');
set(gca, 'xtick', []);
set(gca, 'xticklabel', []);

for k = 1:find((segs(1,:)~=0),1,'last')
    if mod(k,2) % Deal with the case in which the recording starts with a bad segment
        b(k).CData = [0, 0.4470, 0.7410];
    else
        b(k).CData = [1,0,0];
    end
end    
yl = ylim;
ylim([0,yl(2)]);
% xlabel('Recording ID');
xlabel('Length of the recording (minutes)');
legend({'Good','Bad'},'Location','southeast');
box off;
title(['Segmentation of the data after bad segment removal for ' regexprep(bidsID, "_", "\\_")]);


% Calculate bad segment time & percentage for the report
bad_segs_index = ~mod(1:length(segs(1,:)), 2);
bad_segs = segs(1,bad_segs_index);
bs_secs = sum(bad_segs(1,:))*60;
bs_pc = sum(bad_segs(1,:))/sum(segs(1,:))*100;

end

