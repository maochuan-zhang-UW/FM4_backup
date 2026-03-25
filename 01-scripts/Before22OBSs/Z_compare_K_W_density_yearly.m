clear; close all;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/K_W_041525.mat');
data = data_W;

% Convert lat/lon to x/y coordinates
[x, y] = latlon2xy(data(:,2), -data(:,3));
[x_K, y_K] = latlon2xy(data_K(:,2), -data_K(:,3));

% Define spatial constraints
x_limits = [-3 3];   % km
y_limits = [-4.5 4.5]; % km

% Define time periods (one for each year from 2016 to 2025)
years = 2016:2025;
periods = cell(length(years), 1);
period_titles = cell(length(years), 1);
for i = 1:length(years)
    if years(i) == 2025
        periods{i} = {'2025-01-01', '2025-04-15'};
    else
        periods{i} = {sprintf('%d-01-01', years(i)), sprintf('%d-12-31', years(i))};
    end
    period_titles{i} = sprintf('Year: %d', years(i));
end

% Create single figure
figure;
set(gcf, 'Position', [100, 100, 600, 2200]); % Adjusted for 10 rows, slightly wider

% Initialize color limits
c_max = 0;

% First pass to determine global color limits
for p = 1:length(periods)
    % Filter data for current period and depth < 0.1 km
    time_mask = (data(:,1) >= datenum(periods{p}{1})) & (data(:,1) <= datenum(periods{p}{2})) & (data(:,4) < 0.1);
    time_mask_K = (data_K(:,1) >= datenum(periods{p}{1})) & (data_K(:,1) <= datenum(periods{p}{2})) & (data_K(:,4) < 0.1);
    
    for catalog = 1:2
        if catalog == 1
            x_vals = x(time_mask);
            y_vals = y(time_mask);
        else
            x_vals = x_K(time_mask_K);
            y_vals = y_K(time_mask_K);
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

% Set logarithmic color limits (avoid log(0) issues)
log_cmax = log10(max(c_max, 1) + 1);

% Second pass to create plots
for p = 1:length(periods)
    % Filter data for current period and depth < 0.1 km
    time_mask = (data(:,1) >= datenum(periods{p}{1})) & (data(:,1) <= datenum(periods{p}{2})) & (data(:,4) < 0.1);
    time_mask_K = (data_K(:,1) >= datenum(periods{p}{1})) & (data_K(:,1) <= datenum(periods{p}{2})) & (data_K(:,4) < 0.1);
    
    for catalog = 1:2
        subplot( 2, length(years),(p-1)*2 + catalog);
        
        if catalog == 1
            % William's catalog
            x_vals = x(time_mask);
            y_vals = y(time_mask);
            catalog_name = 'William Catalog';
        else
            % Kaiwen's catalog
            x_vals = x_K(time_mask_K);
            y_vals = y_K(time_mask_K);
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
        [X, Y] = meshgrid(x_edges(1:end-1)+0.05, y_edges(1:end-1)+0.05);
        
        % Plot with consistent color scale
        h = pcolor(X, Y, log10(N' + 1));
        set(h, 'EdgeColor', 'none');
        hold on;
        
        % Plot caldera rim
        try
            axial_calderaRim;
            [calderaRim(:,2), calderaRim(:,1)] = latlon2xy(calderaRim(:,2), calderaRim(:,1));
            plot3(calderaRim(:,2), calderaRim(:,1), zeros(size(calderaRim(:,1))), 'white', 'LineWidth', 1);
        catch
            warning('Failed to plot caldera rim for year %d, catalog %d', years(p), catalog);
        end
        
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
        
        % Add text if no data
        if isempty(x_vals)
            text(mean(x_limits), mean(y_limits), 'No Data', ...
                'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 10);
        end
    end
end

% Add shared colorbar
h = colorbar('Position', [0.93 0.05 0.02 0.9]);
h.Label.String = 'log_{10}(Count + 1)';

% Add overall title
sgtitle('Epicenter Density Comparison (0.1 km × 0.1 km, Depth < 0.1 km)', 'FontSize', 16, 'FontWeight', 'bold');