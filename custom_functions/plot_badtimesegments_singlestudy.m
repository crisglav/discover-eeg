function [bs_plot, bs_length, recording_length] = plot_badtimesegments_singlestudy(params,bidsID)
if contains(bidsID,'_')
    x = strsplit(bidsID,'_');
    x = x{1};
else
    x = bidsID;
end
datapath = fullfile(params.PreprocessedDataPath,x,'eeg',[bidsID '_eeg.set']);

% Read header file
hdr = ft_read_header(datapath);
etc = hdr.orig.etc;
events = etc.eventsAfterCRD;

nSamples = length(etc.clean_sample_mask);
nevents = length(events);
segs = zeros(2,max(nevents)*2+1);

mask = etc.clean_sample_mask;
badsegs = reshape(find(diff([false ~mask(:)' false])),2,[]);
s = [1; badsegs(:); nSamples];
s = s(2:end)-s(1:end-1); 
segs(1,1:length(s)) = s;
segs = segs./(hdr.orig.srate*60); % Divide by sampling rate

bs_plot = figure('Units','centimeters','Position', [0 0 18 2],'Visible','Off');
b = barh(1, segs(1,:),'stacked','EdgeColor','none');
set(b,'FaceColor','Flat');
for k = 1:find((segs(1,:)~=0),1,'last')
    if mod(k,2) % Deal with the case in which the recording starts with a bad segment
        b(k).CData = [0, 0.4470, 0.7410];
    else
        b(k).CData = [1,0,0];
    end
end

legend({'Good','Bad'},'Location','northeastoutside');
set(gca,'ytick',[],'yticklabel',{''});
box('off')
xlabel('Time (minutes)');
title('Bad time segments');

% Calculate bad segment time & percentage for the report
bs_length = sum(~etc.clean_sample_mask)/hdr.orig.srate;
recording_length = nSamples/hdr.orig.srate;            

% bad_segs_index = ~mod(1:length(segs(1,:)), 2);
% bad_segs = segs(1,bad_segs_index);
% bs_secs = sum(bad_segs(1,:))*60;
% bs_pc = sum(bad_segs(1,:))/sum(segs(1,:))*100;

end

