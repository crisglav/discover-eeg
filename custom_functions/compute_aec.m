function compute_aec(params,bidsID,freqBand)
% Load EEG data
data = load_preprocessed_data(params,bidsID);

% Band-pass filter the data in the relevant frequency band
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = params.freq_band.(freqBand);
data = ft_preprocessing(cfg, data);

% Load source reconstruction data
load(fullfile(params.source_folder,[bidsID '_source_' freqBand '.mat']),'source');

% Reconstruct the virtual time series (apply spatial filter to sensor level
% data)
cfg  = [];
cfg.pos = source.pos(source.inside,:);
virtChan_data = ft_virtualchannel(cfg,data,source);

% Compute AEC on every trial, then average across trials
tic
connMatrix = aecConnectivity_brainstorm(virtChan_data);
connMatrix = mean(connMatrix,3);
t = toc;
disp([bidsID '_', freqBand, ' computation took ', num2str(t/60), ' minutes'])

save(fullfile(params.connectivity_folder,[bidsID '_aec_' freqBand '.mat']),'connMatrix')

% load(fullfile(params.connectivity_folder,[bidsID '_aec_' freqBand]));
% connMatrix = mean(aec,3);
% connMatrix_Stefan = aec(:,:,1);

%% % Source.avg.mom should be the same as filter_array * avgdata.avg
% cfg = [];
% cfg.covariance = 'no';
% cfg.keeptrials = 'no';
% cfg.removemean = 'no';
% avgdata = ft_timelockanalysis(cfg,avgdata);
% 
% mom = filter_array*avgdata.avg;
% mom_orig = cell2mat(source.avg.mom(source.inside));
%% Fieldtrip version
% cfg = [];
% cfg.trials = 1;
% dataSelect = ft_selectdata(cfg, virtChan_data);
% plot(dataSelect.trial{1}')
% 
% % providing data in hilbert format doesnt work because ft_connectivity
% % analysis explicitely expects dimord = 'rpttap_chan_freq'. 
% cfg = [];
% cfg.method = 'hilbert';
% cfg.foilim = [4 8];
% cfg.output = 'fourier';
% cfg.toi = dataSelect.time{1};
% cfg.keeptrials = 'yes';
% cfg.keeptapers  = 'yes';
% freqSelectHT = ft_freqanalysis(cfg, dataSelect);
% 
% cfg = [];
% cfg.hilbert = 'complex';
% dataSelectHT = ft_preprocessing(cfg, dataSelect);
% 
% %compare the two approaches
% viaFreq = squeeze(freqSelectHT.fourierspctrm(1,:,1,:)); % this only takes into account 4 Hz
% viaPreproc = dataSelectHT.trial{1};
% rel_difference = abs(viaFreq - viaPreproc)./abs(viaFreq);
% figure; plot(rel_difference');
% connMatrix_ft = ft_connectivity_powcorr_ortho(viaFreq);

% diffConMatrices = abs(connMatrix_ft-connMatrix_Stefan);
%%
% % Average over trials
% cfg = [];
% cfg.covariance = 'no';
% cfg.keeptrials = 'no';
% cfg.removemean = 'no';
% avgdata = ft_timelockanalysis(cfg,data);

% % Spectral analysis using morlet wavelets
% cfg = [];
% cfg.method = 'wavelet';
% cfg.foi = 4:1:8;
% cfg.width = 5.83; % 1/2 octave
% cfg.toi = 5;
% cfg.output = 'fourier';
% tfr = ft_freqanalysis(cfg,data);
% tfr.dimord = 'rpttap_chan_freq';
% % pow_corr_ortho using freq datatype
% dat = tfr.fourierspctrm(:,:,1).';
% datout = ft_connectivity_powcorr_ortho(dat, 'refindex','all','tapvec',ones(1,size(dat,2)));

% cfg = [];
% cfg.hilbert = 'complex';
% cfg.keeptrials = 'yes';
% HA = ft_preprocessing(cfg,virtChan_data);

% % pow_corr_ortho using source datatype
% dimord = getdimord(source, 'mom');
% dimtok = tokenize(dimord, '_');
% posdim = find(strcmp(dimtok, '{pos}'));
% posdim = 4; % we concatenate across positions...
% rptdim = find(~cellfun('isempty',strfind(dimtok, 'rpt')));
% rptdim = rptdim-1; % the posdim has to be taken into account...
% dat    = cat(4, data.mom{data.inside});
% dat    = permute(dat,[posdim rptdim setdiff(1:ndims(dat),[posdim rptdim])]);
% 
% datout = ft_connectivity_powcorr_ortho(dat, optarg{:});

% This should work but it doesn't because it is hacked
% % computes AEC on the average over trials
% cfg = [];
% cfg.method = 'powcorr_ortho';
% connMatrix = ft_connectivityanalysis(cfg,source);
end