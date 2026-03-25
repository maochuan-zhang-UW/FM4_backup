clear; close all;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/K_W_041525.mat');
data = data_W;

% Convert lat/lon to x/y coordinates
[x,y] = latlon2xy(data(:,2), -data(:,3));
[x_K,y_K] = latlon2xy(data_K(:,2), -data_K(:,3));
depth = data(:,4);       % Depth values
depth_K = data_K(:,4);   % Depth values for Kaiwen's catalog

% Define spatial constraints
x_limits = [-3 3];   % km
y_limits = [-4.5 4.5]; % km

% Define depth slices (km)
depth_slices = 0:0.3:1.5;
depth_titles = {'0.0-0.3 km', '0.3-0.6 km', '0.6-0.9 km', '0.9-1.2 km', '1.2-1.5 km'};

% Define time periods
periods = {
    {'2022-01-01', '2023-12-31'};
    {'2024-01-01', '2025-12-31'}
};
period_titles = {
    sprintf('Period: %s to %s', datestr(datenum(periods{1}{1}), 'yyyy-mmm-dd'), datestr(datenum(periods{1}{2}), 'yyyy-mmm-dd'));
    sprintf('Period: %s to %s', datestr(datenum(periods{2}{1}), 'yyyy-mmm-dd'), datestr(datenum(periods{2}{2}), 'yyyy-mmm-dd'))
};

% Create figures for each depth slice
for d = 1:length(depth_slices)-1
    depth_min = depth_slices(d);
    depth_max = depth_slices(d+1);
    
    figure;
    set(gcf, 'Position', [100, 100,  572   848]);
    
    % Initialize color limits
    c_max = 0;
    
    % First pass to determine global color limits for this depth slice
    for p = 1:length(periods)
        % Time masks
        time_mask = (data(:,1) >= datenum(periods{p}{1})) & (data(:,1) <= datenum(periods{p}{2}));
        time_mask_K = (data_K(:,1) >= datenum(periods{p}{1})) & (data_K(:,1) <= datenum(periods{p}{2}));
        
        % Depth masks
        depth_mask = (depth >= depth_min) & (depth < depth_max);
        depth_mask_K = (depth_K >= depth_min) & (depth_K < depth_max);
        
        for catalog = 1:2
            if catalog == 1
                x_vals = x(time_mask & depth_mask);
                y_vals = y(time_mask & depth_mask);
            else
                x_vals = x_K(time_mask_K & depth_mask_K);
                y_vals = y_K(time_mask_K & depth_mask_K);
            end
            
            % Apply spatial constraints
            spatial_mask = (x_vals >= x_limits(1)) & (x_vals <= x_limits(2)) & ...
                          (y_vals >= y_limits(1)) & (y_vals <= y_limits(2));
            x_vals = x_vals(spatial_mask);
            y_vals = y_vals(spatial_mask);
            
            % Create histogram
            x_edges = x_limits(1):0.1:x_limits(2);
            y_edges = y_limits(1):0.1:y_limits(2);
            N = histcounts2(x_vals, y_vals, x_edges, y_edges);
            
            % Update global color maximum
            current_max = max(N(:));
            if current_max > c_max
                c_max = current_max;
            end
        end
    end
    
    % Set logarithmic color limits
    log_cmax = log10(c_max + 1);
    
    % Second pass to create plots
    for p = 1:length(periods)
        % Time masks
        time_mask = (data(:,1) >= datenum(periods{p}{1})) & (data(:,1) <= datenum(periods{p}{2}));
        time_mask_K = (data_K(:,1) >= datenum(periods{p}{1})) & (data_K(:,1) <= datenum(periods{p}{2}));
        
        % Depth masks
        depth_mask = (depth >= depth_min) & (depth < depth_max);
        depth_mask_K = (depth_K >= depth_min) & (depth_K < depth_max);
        
        for catalog = 1:2
            subplot(2, 2, (p-1)*2 + catalog);
            
            if catalog == 1
                % William's catalog
                x_vals = x(time_mask & depth_mask);
                y_vals = y(time_mask & depth_mask);
                catalog_name = 'William Catalog';
            else
                % Kaiwen's catalog
                x_vals = x_K(time_mask_K & depth_mask_K);
                y_vals = y_K(time_mask_K & depth_mask_K);
                catalog_name = 'Kaiwen Catalog';
            end
            
            % Apply spatial constraints
            spatial_mask = (x_vals >= x_limits(1)) & (x_vals <= x_limits(2)) & ...
                          (y_vals >= y_limits(1)) & (y_vals <= y_limits(2));
            x_vals = x_vals(spatial_mask);
            y_vals = y_vals(spatial_mask);
            
            % Create 2D density
            x_edges = x_limits(1):0.1:x_limits(2);
            y_edges = y_limits(1):0.1:y_limits(2);
            [N, ~, ~] = histcounts2(x_vals, y_vals, x_edges, y_edges);
            
            % Create grid for pcolor
            [X,Y] = meshgrid(x_edges(1:end-1)+0.05, y_edges(1:end-1)+0.05);
            
            % Plot with consistent color scale
            h = pcolor(X, Y, log10(N' + 1));
            set(h, 'EdgeColor', 'none');
            hold on;
            axial_calderaRim;
[calderaRim(:,2), calderaRim(:,1)] = latlon2xy(calderaRim(:,2), calderaRim(:,1));
plot3(calderaRim(:,2), calderaRim(:,1), zeros(size(calderaRim(:,1))), 'white', 'LineWidth', 1);
            
            % Formatting
            colormap(jet);
            caxis([0 log_cmax]);
            title(sprintf('%s\n%s', catalog_name, period_titles{p}));
            xlabel('X (km)');
            ylabel('Y (km)');
            axis equal tight;
            grid on;
            box on;
            xlim(x_limits);
            ylim(y_limits);
        end
    end
    
    % Add shared colorbar
    h = colorbar('Position', [0.93 0.15 0.02 0.7]);
    h.Label.String = 'log_{10}(Count + 1)';
    
    % Add overall title with depth range
    sgtitle(sprintf('Epicenter Density Comparison (%s)\nDepth Range: %s', ...
           '0.1 km × 0.1 km', depth_titles{d}), ...
           'FontSize', 16, 'FontWeight', 'bold');
    savemyfigureFM4(sprintf('Epicenter Density Comparison (%s)\nDepth Range: %s', ...
           '0.1 km × 0.1 km', depth_titles{d}));
end