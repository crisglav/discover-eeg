function data = load_preprocessed_data(params,bidsID)
    
    % Load one subject data
    x = strsplit(bidsID,'_');
    x = x(1:end-1);
    datapath = fullfile(params.preprocessed_data_path,x{:},'eeg',[bidsID '_eeg.set']);
    
    hdr = ft_read_header(datapath); % all the bids information ins contained in the header of the original file
    cfg = [];
    cfg.dataset = datapath;
    data = ft_preprocessing(cfg);
    
    % Patch to select only EEG electrodes (in case you have specified in
    % the electrodes.tsv type 'EEG' in electrodes that were not EEG
    cfg = [];
    cfg.channel = data.elec.label;
    data = ft_selectdata(cfg,data);
end