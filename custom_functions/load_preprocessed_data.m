function data = load_preprocessed_data(params,bidsID)
    
    % Load data
    x = strsplit(bidsID,'_');
    x = x(1:end-1);
    datapath = fullfile(params.preprocessed_data_path,x{:},'eeg',[bidsID '_eeg.set']);
    
    hdr = ft_read_header(datapath); % all the bids information is contained in the header of the original file
    cfg = [];
    cfg.dataset = datapath;
    data = ft_preprocessing(cfg);
    
    % Take the coordinates from the mni template
    elec_template = ft_read_sens(params.elec_template);
    channels = ismember(elec_template.label,data.label);   
    elec.chanpos = elec_template.chanpos(channels,:);
    elec.chantype = elec_template.chantype(channels,:);
    elec.chanunit = elec_template.chanunit(channels,:);
    elec.elecpos = elec_template.elecpos(channels,:);
    elec.label = elec_template.label(channels,:);
    elec.type = 'custom';
    elec.unit = elec_template.unit;

    % Coordinate system of eelgab has to be defined manually
    data.elec.coordsys = 'ctf';
    data.elec = ft_convert_units(data.elec, 'mm');
    data.elec = ft_convert_coordsys(data.elec, 'acpc');

    % electrode positions] 
    figure;
    ft_plot_sens(data.elec,'label','yes')
    ft_plot_axes(data.elec)
    elec = ft_determine_coordsys(data.elec);
    
    % Overlay electrodes with the head model
    load(params.volpath,'vol');
    figure;
    ft_plot_headmodel(vol,'facealpha',0.1,'facecolor',[0.1 0.1 0.1],'edgecolor',[1 1 1],'edgealpha',0.5);
    hold on;
    ft_plot_sens(elec,'style','r','label','label','elec','true','elecshape','disc','elecsize',5,'facecolor','r');
    view(90,0);


    % Patch to select only EEG electrodes (in case you have specified in
    % the electrodes.tsv type 'EEG' in electrodes that were not EEG
    cfg = [];
    cfg.channel = data.elec.label;
    data = ft_selectdata(cfg,data);
end