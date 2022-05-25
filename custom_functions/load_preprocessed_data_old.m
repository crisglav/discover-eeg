function data = load_preprocessed_data_old(params,bidsID)
    
    % Load one subject data
    s = strsplit(bidsID,'_');
    datapath = fullfile(params.bids_folder,s{1},'eeg',[bidsID '_eeg.set']);
    
    hdr = ft_read_header(datapath);
    cfg = [];
    cfg.dataset = datapath;
    data = ft_preprocessing(cfg);
    
    % Segment the data into epochs
    cfg = [];
    cfg.length = params.epoch_length;
    cfg.overlap = params.epoch_overlap;
    data = ft_redefinetrial(cfg,data);
    
    % Reject epochs that contain discontinuities
    discontinuities = hdr.orig.event(strcmp('boundary', {hdr.orig.event.type}));
    clear hdr;
    if isempty(discontinuities)
        badIntervalMatrix = [];
    else
        badIntervalMatrix(:,1) = [discontinuities.latency];
        badIntervalMatrix(:,2) = [discontinuities.latency];
    end
    cfg = [];
    cfg.artfctdef.reject  = 'complete';
    cfg.artfctdef.visual.artifact = badIntervalMatrix;
    data = ft_rejectartifact(cfg,data);
    
%     % Electrode layout (only clean channels)
%     clean_channels = ismember(params.elec_template.label,data.label);
%     elec.chanpos = elec_template.chanpos(clean_channels,:);
%     elec.chantype = elec_template.chantype(clean_channels,:);
%     elec.chanunit = elec_template.chanunit(clean_channels,:);
%     elec.elecpos = elec_template.elecpos(clean_channels,:);
%     elec.label = elec_template.label(clean_channels,:);
%     elec.type = 'custom';
%     elec.unit = elec_template.unit;
end