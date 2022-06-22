function data = load_preprocessed_data(params,bidsID)
    
    % Load data
    x = strsplit(bidsID,'_');
    x = x(1:end-1);
    datapath = fullfile(params.preprocessed_data_path,x{:},'eeg',[bidsID '_eeg.set']);
    
    hdr = ft_read_header(datapath); % all the bids information is contained in the header of the original file
    cfg = [];
    cfg.dataset = datapath;
    data = ft_preprocessing(cfg);
    
    % Using electrode positions from EEGLab
    if strcmp(params.bidschanloc,'on')
%         % Coordinate system of eelgab has to be defined manually
%         data.elec.coordsys = 'ctf';
%         data.elec = ft_convert_units(data.elec, 'mm');
%         data.elec = ft_convert_coordsys(data.elec, 'acpc');
%         
%         % electrode positions
%         figure;
%         ft_plot_sens(data.elec,'label','yes')
%         ft_plot_axes(data.elec)
%         
%         % Overlay electrode positions from EEGLab with the head model -> not
%         % properly aligned. It is better to define the electrode positions
%         % directly from the aligned template (see below)
%         load(params.volpath,'vol');
%         figure;
%         ft_plot_headmodel(vol,'facealpha',0.1,'facecolor',[0.1 0.1 0.1],'edgecolor',[1 1 1],'edgealpha',0.5);
%         hold on;
%         ft_plot_sens(data.elec,'style','r','label','label','elec','true','elecshape','disc','elecsize',5,'facecolor','r');
%         view(90,0);
        
        elecpath = fullfile(params.raw_data_path,x{:},'eeg',[x{:} '_electrodes.tsv']);
        elec_tsv = ft_read_sens(elecpath);
        channels = ismember(elec_tsv.label,data.label);
               
    else    
        % Take the coordinates from the mni template
        elec_template = ft_read_sens('standard_1005.elc');
        channels = ismember(elec_template.label,data.label);   
        elec.chanpos = elec_template.chanpos(channels,:);
        elec.chantype = elec_template.chantype(channels,:);
        elec.chanunit = elec_template.chanunit(channels,:);
        elec.elecpos = elec_template.elecpos(channels,:);
        elec.label = elec_template.label(channels,:);
        
        % Arrange the electrodes as in the original data
        c = setdiff(data.label,elec_template.label);
        ic = ismember(data.label,elec_template.label);
        c = sprintf('%s ', c{:});
        warning(['The position of the channel(s) ' c 'was not found in a standard template and they will be removed. If you want to include them please specify their position in a electrodes.tsv and change define_params accordingly.']);           
        
        cfg = [];
        cfg.channel = data.elec.label(ic);
        data = ft_selectdata(cfg,data);
        % TO DO:
        [a,i] = sort(data.label);
        elec.label(i);
        
    
    end
    elec.chanpos = elec_template.chanpos(channels,:);
    elec.chantype = elec_template.chantype(channels,:);
    elec.chanunit = elec_template.chanunit(channels,:);
    elec.elecpos = elec_template.elecpos(channels,:);
    elec.label = elec_template.label(channels,:);
    elec.type = 'custom';
    elec.unit = elec_template.unit;
    
    % Overlay electrode positions
    load(params.volpath,'vol');
    figure;
    ft_plot_headmodel(vol,'facealpha',0.1,'facecolor',[0.1 0.1 0.1],'edgecolor',[1 1 1],'edgealpha',0.5);
    hold on;
    ft_plot_sens(elec_tsv,'label','number','elec','true','elecshape','disc','elecsize',5,'facecolor','r');
    view(90,0);
    
    data.elec = elec;
    
%     % Patch to select only EEG electrodes (in case you have specified in
%     % the electrodes.tsv type 'EEG' in electrodes that were not EEG)
%     cfg = [];
%     cfg.channel = data.elec.label;
%     data = ft_selectdata(cfg,data);
end