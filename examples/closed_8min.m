function [trl, event] = closed_8min(cfg)
% Select eyes closed segments (marker 'S210') 

hdr = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);
samples = {event.sample};

% Detect when the s210 segments start and finish
s210 = diff([0 strcmp({event.value},'S210') 0]);
% Check that there are S210 markers
assert(any(s210),'No s210 marker was found.');

start_samples = [samples{s210 == 1}];
end_samples = [samples{find(s210 == -1)-1}];

% % Find the segments longer than 50s
% ix = find((end_samples-start_samples),1);
%  % Check that there are at least a segment longer than 50s
% assert(~isempty(ix),'No segment longer than 50s found.');

trl(:,1) = start_samples;
trl(:,2) = end_samples;
trl(:,3) = 0;       
   
end
