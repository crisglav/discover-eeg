%% rm_raincloud_cg - plots a raincloud of different groups and observations
% This function is based on the rm_raincloud function.m of
% github.com/RainCloudPlots/RaincloudPlots/
%
% Please see the publication Allen M, Poggliali D, Whitaker K et al.
% Raincloud plots: a multi-platform tool for robust data visualization
% [version 2; peer revew:2 approved]. Wellcome Open Res 2021, 4:63. DOI:
% 10.12688/wellcomeopenres.15191.2
%
% Use like: h = rm_raincloud_cg(data,'parameter','value')
% Where 'data' is an M x N cell array of M measurements and N data series
% See below for optional parameters.
% 
% ----------------------------- NEW! ------------------------------------
% - Modfified 'colours' parameter. Now an array of colours can be given so
% that each measure has a different colour.
% - Renamed 'alpha' parameter to 'opacity'. Controls the opacity of the
% cloud and the raindrops. Range 0 - 1
% - Added 'box_width' parameter. Sets the width of the boxplots and rain
% jitter to a value specified by the user.
% - Deleted 'plot_top_to_bottom' parameter
% - Changed logical parameters: now input must be true or false (before 0 or 1)
% - Added limits to y axis so that plots remain aligned
% - Fixed typos: mean to median
%
% Modified: Cristina Gil, 26.03.2020
% 
% - Fixed computation of boxplots
% - Default colors are now matlab standard colors
% Modified: Cristina Gil, 27.09.2021

function h = rm_raincloud_cg(data, varargin)
%% ---------------------------- INPUT ----------------------------
%
% data - M x N cell array of M measurements and N data series
%
% --------------------- OPTIONAL ARGUMENTS ----------------------
%
% colours               - N x 3 x M array defining the colour to plot each data series
% density_type          - choice of density algo ('ks' or 'rath'). Default = 'ks'
% bandwidth             - Controls the smooth of the cloud. Set to empty [] to automatically adjust to the data (default = 1)
% plot_median_lines     - logical. Set to false if you don't want median lines plotted (default = true)
% plot_median_dots      - logical. Set to false if you don't want median dots plotted (default = false)
% box_on                - logical. Turns boxplots on/off (default = true)
% line_width            - scalar value to set global line width (default = 2)
% box_width             - scalar value that defines the width of the boxplot and the rain (default = 0.1)
% bxcl                  - Array double 3x1. Color of box outline (defalut [0 0 0])
% box_col_match         - logical. If true boxes match the colour of clouds (default = false)
% box_dodge             - logical. Turn on/off box plot dodging (default = false)
% raindrop_size         - scalar. Positive value to control the size of the raindrops (default = 3)
% opacity               - double with range 0 to 1. Controls the opacity of the cloud (defalut = 0.5)
% dist_plots            - scalar. Defines the distance between plots (default = 1.5)  
% aligned_plots         - logical. Align plots in the same x axes (default = true). 

% ---------------------------- OUTPUT ----------------------------
% h is a cell array containing handles of the various figure parts:
% h.p{i,j}  is the handle to the density plot from data{i,j}
% h.s{i,j}  is the handle to the 'raindrops' (individual datapoints) from data{i,j}
% h.b1{i,j} [optional: only if box_on is true] is the handle for the mean line of the boxplots
% h.b2{i,j} [optional: only if box_on is true] is the handle for right whisker of the boxplots
% h.b3{i,j} [optional: only if box_on is true] is the handle for left whisker of the boxplots
% h.m(i,j)  [optional: only if plot_median_dots is true] is the handle to the single, large dot that represents median(data{i,j})
% h.l(i,j)  [optional: only if plot_median_lines is true] is the handle for the line connecting h.m(i,j) and h.m(i+1,j)
%
% ------------------------ EXAMPLE USAGE -------------------------
% h = rm_raincloud_cg(raincloudData,'colours',[0 .3961 .7412; .8902 .4471 .1333],'plot_mean_lines',0,'box_on',1)
% 
%% check dimensions of data
% input parsing settings
p = inputParser;
p.CaseSensitive = true;
p.Parameters;
p.Results;
p.KeepUnmatched = true;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);

[n_plots_per_series, n_series] = size(data);
cmap = colormap(lines);
col = cmap(1:n_series,:);
col = repmat(col,[1,1,n_plots_per_series]);
col = permute(col,[3 2 1]);
%% default arguments
% set the desired and optional input arguments
addRequired(p, 'data', @iscell);
addOptional(p, 'colours', col, @isnumeric)
addOptional(p, 'density_type', 'ks', @ischar)
addOptional(p, 'bandwidth', 1)
addOptional(p, 'plot_median_lines', true, @islogical)
addOptional(p, 'plot_median_dots', false, @islogical)
addOptional(p, 'box_on', true, @islogical)
addOptional(p, 'bxcl', [0 0 0], @isnumeric)
addOptional(p, 'line_width', 2, validScalarPosNum)
addOptional(p, 'box_width',0.1, @isnumeric)
addOptional(p, 'box_col_match', false, @islogical)
addOptional(p, 'box_dodge', false, @islogical)
addOptional(p, 'raindrop_size', 3, validScalarPosNum)
addOptional(p, 'opacity', 0.5, @isnumeric)
addOptional(p, 'dist_plots', 1.5, @isnumeric)
addOptional(p, 'aligned_plots',true, @islogical)


% parse the input
parse(p,data,varargin{:});
% then set/get all the inputs out of this structure
data                = p.Results.data;
colours             = p.Results.colours;
bandwidth           = p.Results.bandwidth;
density_type        = p.Results.density_type;
plot_median_lines   = p.Results.plot_median_lines;
plot_median_dots    = p.Results.plot_median_dots;
box_on              = p.Results.box_on;
bxcl                = p.Results.bxcl;
line_width          = p.Results.line_width;
box_width           = p.Results.box_width;
box_col_match       = p.Results.box_col_match;
box_dodge           = p.Results.box_dodge;
raindrop_size       = p.Results.raindrop_size;
opacity             = p.Results.opacity;
dist_plots          = p.Results.dist_plots;
aligned_plots       = p.Results.aligned_plots;

%% Calculate properties of density plots

% Probably okay to hard-code this as it just determines the granularity of
% the density estimate
density_granularity = 200;

n_bins = repmat(density_granularity, n_plots_per_series, n_series);

% initialize variables
ks = cell(size(data));
x = cell(size(data));
q = cell(size(data));
faces = cell(size(data));

% calculate kernel densities
for i = 1:n_plots_per_series
    for j = 1:n_series
       
        switch density_type
            
            case 'ks'
                
                % compute density using 'ksdensity'
                [ks{i, j}, x{i, j}] = ksdensity(data{i, j}, 'NumPoints', n_bins(i, j), 'bandwidth', bandwidth, 'Function','pdf');
                
            case 'rash'
                
                % check for rst_RASH function (from Robust stats toolbox) in path, fail if not found 
                assert(exist('rst_RASH', 'file') == 2, 'Could not compute density using RASH method. Do you have the Robust Stats toolbox on your path?');
                
                % compute density using RASH
                [x{i, j}, ks{i, j}] = rst_RASH(data{i, j});
                
                % override default 'n_bins' as rst_RASH determines number of bins
                n_bins(i, j) = size(ks{i, j}, 2);
        end
        
        % Define the faces to connect each adjacent f(x) and the corresponding points at y = 0.
        q{i, j}     = (1:n_bins(i, j) - 1)';
        faces{i, j} = [q{i, j}, q{i, j} + 1, q{i, j} + n_bins(i, j) + 1, q{i, j} + n_bins(i, j)];
        
    end
end

if aligned_plots==1
    spacing = 1/dist_plots;
    ks_offsets = (0:n_plots_per_series-1).*spacing; % For aligned version
    jit_width = box_width*spacing;
else
    % determine spacing between plots
    plotting_space = mean(mean(cellfun(@max, ks)));
    jit_width = plotting_space/8;
    spacing = plotting_space * (n_series+dist_plots); % dist_plots to have a margin. set higher if you want more distance between plots
    ks_offsets = (0:n_plots_per_series-1) .* spacing;
end

% flip so first plot in series is plotted on the *top*
ks_offsets  = fliplr(ks_offsets);

% calculate patch vertices from kernel density
verts = cell(size(data));
for i = 1:n_plots_per_series
    for j = 1:n_series
        verts{i, j} = [x{i, j}', ks{i, j}' + ks_offsets(i); x{i, j}', ones(n_bins(i, j), 1) * ks_offsets(i)];
    end
end

%% boxplot parameters 
if box_on
    Y = cell(size(data));
    for i = 1:n_plots_per_series
        for j = 1:n_series
            quartiles   = quantile(data{i,j},[0.25 0.75 0.5]);
            iqr         = quartiles(2) - quartiles(1);
            Xs          = sort(data{i,j});
%             whiskers(1) = min(Xs(Xs > (quartiles(1) - (1.5 * iqr))));
%             whiskers(2) = max(Xs(Xs < (quartiles(2) + (1.5 * iqr))));
            % Option B
            if (quartiles(1) - 1.5 * iqr > min(Xs))
                whiskers(1) = min(Xs(Xs > (quartiles(1) - (1.5 * iqr))));
            else
                whiskers(1) = min(Xs);
            end
            
            if (quartiles(2) + 1.5 * iqr < max(Xs))
                whiskers(2) = max(Xs(Xs < (quartiles(2) + (1.5 * iqr))));
            else
                whiskers (2) = max(Xs);
            end
            Y{i,j}      = [quartiles whiskers];
           
        end
    end
end

%% jitter for the raindrops
drop_pos = cell(size(data));
for i = 1:n_plots_per_series
    for j = 1:n_series
        if box_dodge
            jit = rand(1, length(data{i,j}))*jit_width;
            offset = ks_offsets(i)-(4*j-1)*jit_width;
            drop_pos{i,j}=offset-jit;
        else
            jit = rand(1, length(data{i,j}))*jit_width;
            offset = ks_offsets(i)-(4*j-2.5)*jit_width;
            drop_pos{i,j}=offset-jit;

        end
    end
end

%% plot cloud, raindrop and boxplot
hold on

% patches
offsets = zeros(n_plots_per_series,n_series);
for i = 1:n_plots_per_series
    for j = 1:n_series
        
        % plot patches
        h.p{i, j} = patch('Faces', faces{i, j}, 'Vertices', verts{i, j}, 'FaceVertexCData', colours(j, :, i), 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', opacity);
        
        % scatter rainclouds
        h.s{i, j} = scatter(data{i, j}, drop_pos{i,j}, 'MarkerFaceColor', colours(j, :, i), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', opacity, 'SizeData', raindrop_size);
    
        % Plot boxplots
        if box_on
            if box_col_match
                bxcl = colours(j,:,i);
            end
            if box_dodge
                offsets(i,j) = ks_offsets(i)-(4*j-2.5)*jit_width;
            else
                offsets(i,j) = ks_offsets(i)-(4*j-2)*jit_width;
            end
                box_pos = [Y{i,j}(1) offsets(i,j)-jit_width*0.5 Y{i,j}(2)-Y{i,j}(1) jit_width];
                % median line
                h.b1{i,j} = line([Y{i,j}(3) Y{i,j}(3)], [offsets(i,j)-jit_width*0.5 offsets(i,j)+jit_width*0.5], 'col', bxcl, 'LineWidth', line_width);
                % whiskers
                h.b2{i,j} = line([Y{i,j}(2) Y{i,j}(5)], [offsets(i,j) offsets(i,j)], 'col', bxcl, 'LineWidth', line_width);
                h.b3{i,j} = line([Y{i,j}(1) Y{i,j}(4)], [offsets(i,j) offsets(i,j)], 'col', bxcl, 'LineWidth', line_width);
                
                % 'box' of the 'boxplot'
                h.b4{i,j} = rectangle('Position', box_pos);
                set(h.b4{i,j}, 'EdgeColor', bxcl)
                set(h.b4{i,j}, 'LineWidth', line_width);
        
        end
    end
end

%% plot median dots/lines
cell_medians = cellfun(@median, data);

% plot median lines
if plot_median_lines
    for i = 1:n_plots_per_series - 1 % We have n_plots_per_series-1 lines because lines connect pairs of points
        for j = 1:n_series
            if box_col_match
                bxcl = colours(j,:,i);
            end
            % median line reaches the boxplot
            h.l(i, j) = line(cell_medians([i i+1], j), [offsets(i,j) offsets(i+1,j)], 'LineWidth', line_width, 'Color', bxcl);
        end
    end
end

% plot median dots
if plot_median_dots
    for i = 1:n_plots_per_series
        for j = 1:n_series
            h.m(i, j) = scatter(cell_medians(i, j), ks_offsets(i), 'MarkerFaceColor', colours(j, :, i), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 1, 'SizeData', raindrop_size * 5, 'LineWidth', 2);
        end
    end
end

%% clear up axis labels
% 'YTick', likes values that *increase* as you go up the Y-axis, but we plot the first
% raincloud at the top. So flip the vector around
set(gca, 'YTick', fliplr(ks_offsets));
ylim([-spacing/1.5 , spacing*n_plots_per_series-spacing+spacing/1.5])
set(gca, 'YTickLabel', n_plots_per_series:-1:1);

%% plot rotation
% NOTE: Because it's easier, we actually prepare everything plotted
% top-to-bottom, then - by default - we rotate it here.
view([90 -90]);
axis ij
