% EEGLAB history file generated on the 15-Jun-2022
% ------------------------------------------------
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
[EEG ALLEEG CURRENTSET] = eeg_retrieve(ALLEEG,1);
EEG = eeg_checkset( EEG );
figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
EEG = eeg_checkset( EEG );
figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab('rebuild');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
eeglab('redraw');
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
EEG = eeg_checkset( EEG );
EEG=pop_chanedit(EEG, 'lookup','/rechenmagd4/toolboxes_and_functions/eeglab/plugins/dipfit4.3/standard_BEM/elec/standard_1005.elc','changefield',{31,'X','-65'},'changefield',{31,'Y','50'},'changefield',{31,'Z','-60'},'changefield',{32,'X','65'},'changefield',{32,'Y','50'},'changefield',{32,'Z','-60'});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
pop_saveh( EEG.history, 'eeglabhist.m', '/rechenmagd4/Experiments/2021_preprocessing/pipeline/');
eeglab redraw;
