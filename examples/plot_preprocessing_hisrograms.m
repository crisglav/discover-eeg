% Load preprocessed data
params = define_params();
[STUDY, ~] = pop_loadstudy('filename', [params.study '-clean.study'], 'filepath', params.preprocessed_data_path);

etc = {ALLEEG.etc};
nRec = length(ALLEEG);

% Bad channels histogram
recmask = cellfun(@(x) isfield(x, 'clean_channel_mask'), etc); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_channel_mask, etc(recmask),'uni',0);
tmp = cellfun(@(x) x(:), tmp, 'UniformOutput', false); % force to be a column vector
badrecs = double(~cat(2,tmp{:}));
badchans = sum(badrecs,1);
f = figure;
histogram(badchans);
ylabel('number of recordings');
xlabel('number of bad channels');
saveas(f,fullfile(params.preprocessed_data_path,'preprocessing_visualization','badchan_hist.jpg'));
mean(badchans)
std(badchans)

% Bad ICs histogram
extract_bad_ICs = @(x) sum((x(:,2:end) > params.IClabel(2:end,1)').*(x(:,2:end) < params.IClabel(2:end,2)'));
classifications = cellfun(@(x) x.ic_classification.ICLabel.orig_classifications, etc,'uni',0);
badics = cellfun(extract_bad_ICs, classifications,'uni',0);
classes = etc{1}.ic_classification.ICLabel.classes;
badics = reshape(cell2mat(badics),length(classes)-1,nRec);
badics = sum(badics,1);
f = figure;
histogram(badics);
ylabel('number of recordings');
xlabel('number of bad ICs');
saveas(f,fullfile(params.preprocessed_data_path,'preprocessing_visualization','badICs_hist.jpg'));
mean(badics)
std(badics)

% Bad segments histogram
srates = cell2mat({ALLEEG.srate})';
recmask = cellfun(@(x) isfield(x, 'clean_sample_mask'), etc); % Deal with the case where no bad channels were detected
tmp = cellfun(@(x) x.clean_sample_mask, etc(recmask),'uni',0);
badsegs = cellfun(@(x) sum(~x), tmp)./(srates');
f = figure;
histogram(badsegs,'BinEdges',linspace(0,200,11));
ylabel('number of recordings');
xlabel('total amount of rejected segments (seconds)');
saveas(f,fullfile(params.preprocessed_data_path,'preprocessing_visualization','badsegs_hist.jpg'));
mean(badsegs)
std(badsegs)

% Clean epochs left histogram
nepochs = [ALLEEG.trials];
f = figure;
histogram(nepochs,'BinEdges',0:10:80);
ylabel('number of recordings');
xlabel('number of clean epochs');
saveas(f,fullfile(params.preprocessed_data_path,'preprocessing_visualization','epochs_hist.jpg'));
mean(nepochs)
std(nepochs)