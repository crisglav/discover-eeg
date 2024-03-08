function data = load_preprocessed_data(params,bidsID)
if exist(fullfile(params.PreprocessedDataPath,bidsID,'eeg',[bidsID '_eeg.mat']),'file')
    load(fullfile(params.PreprocessedDataPath,bidsID,'eeg',[bidsID '_eeg.mat']),'data');
else

    % Load data
    ses ='';
    if contains(bidsID,'_')
        bidsparts = strsplit(bidsID,'_');
        sub = bidsparts{1};
        ses_mask = cellfun(@(x) contains(x, 'ses-'), bidsparts);
        if any(ses_mask)
            ses = bidsparts{ses_mask};
        end
    else
        sub = bidsID;
    end
    
    datapath = fullfile(params.PreprocessedDataPath,sub,ses,'eeg',[bidsID '_eeg.set']);
    
    hdr = ft_read_header(datapath); % all the bids information is contained in the header of the original file
    cfg = [];
    cfg.dataset = datapath;
    data = ft_preprocessing(cfg);
    
    % Using electrode positions from electrodes.tsv
    if strcmp(params.BidsChanloc,'on')
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
%         load(params.HeadModelPath,'vol');
%         figure;
%         ft_plot_headmodel(vol,'facealpha',0.1,'facecolor',[0.1 0.1 0.1],'edgecolor',[1 1 1],'edgealpha',0.5);
%         hold on;
%         ft_plot_sens(data.elec,'style','r','label','label','elec','true','elecshape','disc','elecsize',5,'facecolor','r');
%         view(90,0);
        if ~isempty(ses)
            elec = ft_read_sens(fullfile(params.RawDataPath,sub,ses,'eeg',[sub '_' ses '_electrodes.tsv']));
        else
            elec = ft_read_sens(fullfile(params.RawDataPath,sub,'eeg',[sub '_electrodes.tsv']));
        end
        channels = ismember(elec.label,data.label);
        
        c = setdiff(elec.label,data.label);
        if ~isempty(c)
            c = sprintf('%s ', c{:});
            warning(['Channels(s) ' c 'are not present in the data and will be discarded']);
        end
        
        data.elec.chanpos(channels,:) = elec.chanpos(channels,:);
        data.elec.chantype(channels,:) = elec.chantype(channels,:);
        data.elec.chanunit(channels,:) = elec.chanunit(channels,:);
        data.elec.elecpos(channels,:) = elec.elecpos(channels,:);
        data.elec.label(channels,:) = elec.label(channels,:);
        data.elec.type = 'custom';
        data.elec.unit = elec.unit;
        
        
        % Check whether there are channels not present in the tsv file
        c = setdiff(data.label,elec.label);
        if ~isempty(c)
            warning(['The position of the channel(s) ' sprintf('%s ', c{:}) 'was not found in ' sub '_electrodes.tsv']);            
        end
        if ismember(c,hdr.orig.ref)
            warning(['Adding the position of the reference ' hdr.orig.ref ' manually, based on define_params.m']);
            ref = ismember(data.label,hdr.orig.ref);
            data.elec.chanpos(ref,:) = [params.RefCoord.X, params.RefCoord.Y, params.RefCoord.Z];
            data.elec.elecpos(ref,:) = [params.RefCoord.X, params.RefCoord.Y, params.RefCoord.Z];
        end
        
    else    
        % EEGLab uses their own coordinate system. When we looked for the
        % electrode positions in EEGLab from the MNI template, then they
        % were transformed EEGlab coordinate system. That is the reason why
        % the values of data.elec differ from elec_template.
        % Best practice is to take the coordinates of the original
        % template.
        
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
        
        if isequal(data.label,elec.label)
            data.elec = elec;
        else
            warning('data.label is not equal to elec.label')
        end
        
    end

    save(fullfile(params.PreprocessedDataPath,sub,ses,'eeg',[bidsID '_eeg.mat']));
end
end