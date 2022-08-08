function data = load_preprocessed_data(params,bidsID)
    
    % Load data
    if contains(bidsID,'_')
        x = strsplit(bidsID,'_');
        x = x{1};
    else
        x = bidsID;
    end
    datapath = fullfile(params.preprocessed_data_path,x,'eeg',[bidsID '_eeg.set']);
    
    hdr = ft_read_header(datapath); % all the bids information is contained in the header of the original file
    cfg = [];
    cfg.dataset = datapath;
    data = ft_preprocessing(cfg);
    
    % Using electrode positions from electrodes.tsv
    if strcmp(params.bidschanloc,'on')
          % Coordinates taken from the .set file are not aligned with mni template!! 
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
        
        elec = ft_read_sens(fullfile(params.raw_data_path,x{:},'eeg',[x{:} '_electrodes.tsv']));
        channels = ismember(elec.label,data.label);
        
        % Check whether there are channels not present in the tsv file
        c = setdiff(data.label,elec.label);
        if ~isempty(c)
            c = sprintf('%s ', c{:});
            warning(['The position of the channel(s) ' c ' was not found in ' x{:} '_electrodes.tsv']);
        end
        
        c = setdiff(elec.label,data.label);
        if ~isempty(c)
            c = sprintf('%s ', c{:});
            warning(['Channels(s) ' c ' are not present in the data and will be discarded']);
        end
        
        elec.chanpos = elec.chanpos(channels,:);
        elec.chantype = elec.chantype(channels,:);
        elec.chanunit = elec.chanunit(channels,:);
        elec.elecpos = elec.elecpos(channels,:);
        elec.label = elec.label(channels,:);
        elec.type = 'custom';
        elec.unit = elec.unit;
        
    else    
        % Take the coordinates from the mni template
        elec_template = ft_read_sens('standard_1005.elc');
        channels = ismember(elec_template.label,data.label);   
        
        chanpos = elec_template.chanpos(channels,:);
        chantype = elec_template.chantype(channels,:);
        chanunit = elec_template.chanunit(channels,:);
        elecpos = elec_template.elecpos(channels,:);
        label = elec_template.label(channels,:);
        
        % Check whether there are non-standard channels
        c = setdiff(data.label,elec_template.label);
        if ~isempty(c)
            c = sprintf('%s ', c{:});
            warning(['The position of the channel(s) ' c 'was not found in a standard template and they will be removed.']);
        end                
        % Remove non-standard channels
        ic = ismember(data.label,elec_template.label);       
        if ~all(ic)
            cfg = [];
            cfg.channel = data.elec.label(ic);
            data = ft_selectdata(cfg,data);
        end
        
        % Arrange the electrodes as in the original data
        [~,ix] = sort(data.label);
        [~,jx] = sort(ix);        
        [~,kx] = sort(label);
        
        elec.chanpos = chanpos(kx(jx),:);
        elec.chantype = chantype(kx(jx),:);
        elec.chanunit = chanunit(kx(jx),:);
        elec.elecpos = elecpos(kx(jx),:);
        elec.label = label(kx(jx),:);
        elec.type = 'custom';
        elec.unit = elec_template.unit;
        
    end
    
    if isequal(data.label,elec.label)
        data.elec = elec;
    else
        warning('data.label is not equal to elec.label')
    end
    
end