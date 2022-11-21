%% An example script to visualize the raw data timeseries programatically in EEGlab
%
% Cristina Gil, 21.11.2022, Technical University of Munich

% Add EEGLab
eeglab_path = '/rechenmagd4/toolboxes_and_functions/eeglab';
run(fullfile(eeglab_path,'eeglab.m'));

%% Plot raw data from single recording
% This would be similar for preprocessed data too.
orig_folder = '/rechenmagd4/Cristina/cpCGX_BIDS/derivatives_v2022_11_15_orig/';
subjID = 'sub-043';
task = 'EC';
filename = [subjID '_task-' task '_eeg.set'];
filepath = fullfile(orig_folder,subjID,'eeg');
EEGRaw = pop_loadset('filename',filename,'filepath',filepath);
% High pass filter the data at 1 Hz for visulization
[EEGRaw,~,~] = pop_eegfiltnew(EEGRaw, 1, []);
pop_eegplot( EEGRaw, 1, 1, 1);


%% Plot raw data timeseries (Before BIDS conversion)
filepath = '/rechenmagd4/Experiments/2022_chronic_pain_CGX/data';
files = dir(fullfile(filepath,'*.vhdr'));
filename = '059ZIS20220922EC.vhdr';
EEGRaw = pop_loadbv(filepath, filename);
% High pass filter the data at 1 Hz for visulization
[EEGRaw,~,~] = pop_eegfiltnew(EEGRaw, 1, []);
pop_eegplot( EEGRaw, 1, 1, 1);