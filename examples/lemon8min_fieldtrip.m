%% Extract 8 minute of resting-state eyes closed from the LEMON dataset
%
% This script extracts eight minutes of eyes-closed resting state from
% the LEMON dataset(concatenated blocks of 1 min duration). You need to have fieldtrip.
%
% Cristina Gil, 11.08.2022, Technical University of Munich

clear all;

% Add fieldtrip
addpath /rechenmagd4/toolboxes_and_functions/fieldtrip
ft_defaults
    
% Specify the folders with the original data and the output directory
dir_in = '/rechenmagd4/Experiments/2021_preprocessing/datasets/LEMON/EEG_Raw_BIDS_ID';
dir_out = '/rechenmagd4/Experiments/2021_preprocessing/datasets/LEMON-8min';
files = dir(fullfile(dir_in,'sub-*'));

for iRec=1:length(files)
    datafile = dir(fullfile(files(iRec).folder,files(iRec).name,'eeg','*.vhdr'));
    
    try
    cfg = [];
    cfg.dataset = fullfile(datafile.folder,datafile.name);
    cfg.trialdef.eventtype = 'Stimulus';
    cfg.trialdef.eventvalue = 'S210';
    cfg.trialfun = 'closed_8min';
    cfg = ft_definetrial(cfg);
    
    data_seg = ft_preprocessing(cfg);
        
    % Concatenate trials
    data_cat = data_seg;
    data_cat.trial = [data_seg.trial{:}];
    data_cat.time = (0:size(data_cat.trial,2)-1)/data_cat.fsample;
    
    % Edit header and event files
    hdr = ft_read_header(cfg.dataset);
    hdr.nSamples = size(data_cat.trial,2);
    
    segment_length = cellfun(@length, data_seg.time);
    new_markers = [1, cumsum(segment_length(1:end-1))];
    event = struct('type',repmat({'S210'},1,length(new_markers)),...
                  'value',repmat({'boundary'},1,length(new_markers)),...
        'sample',num2cell(new_markers),'duration',num2cell(ones(1,length(new_markers))));
    
%     % Plot original data
%     cfg = [];
%     cfg.dataset = fullfile(datafile.folder,datafile.name);
%     data_orig = ft_preprocessing(cfg);
%     figure;
%     plot(data_orig.time{1,1},data_orig.trial{1,1}(1,:));
%     for i=1:8
%         hold on
%         xline(data_seg.sampleinfo(i,1)/2500,'g')
%         hold on
%         xline(data_seg.sampleinfo(i,2)/2500,'m')
%     end
%     
%     % Plot only the eyes closed data with concatenated blocks
%     figure;
%     plot(data_cat.time,data_cat.trial(1,:));
%     for i=1:8
%     hold on
%     xline(new_markers(i)/2500,'g');
%     end
    
    % Write file
    datafileout = fullfile(dir_out,files(iRec).name,'eeg',[files(iRec).name '_eeg.eeg']);
    ft_write_data(datafileout,data_cat.trial,'header',hdr,'event',event,'dataformat','brainvision_eeg');
    catch ME
        warning(['Problem with ' files(iRec).name ': ' ME.message ' Skipping recording'])
    end
end


files_out = dir(fullfile(dir_out,'sub-*'));
sub_in = {files.name};
sub_out = {files_out.name};

sub_in(~contains(sub_in,sub_out))