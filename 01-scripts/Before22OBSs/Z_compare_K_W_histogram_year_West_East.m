clear; close all;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/K_W_041525.mat');
data = data_W;

% Convert lat/lon to x/y coordinates for both datasets
[x,y] = latlon2xy(data(:,2), -data(:,3));
[x_K,y_K] = latlon2xy(data_K(:,2), -data_K(:,3));

% Define parameters
bin_edges = 0:0.1:2;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
years = 2015:2025;
num_years = length(years);

% Create two figures - one for each region
for region = 1:2
    figure;
    set(gcf, 'Position', [100+(region-1)*600, 100, 1200, 900]);
    
    % Determine subplot arrangement (3x4 grid for 11 years)
    rows = 3;
    cols = 4;
    
    for y = 1:num_years
        subplot(rows, cols, y);
        
        % Get data for current year
        start_date = datenum(years(y), 1, 1);
        end_date = datenum(years(y), 12, 31);
        
        % Filter both datasets based on region
        if region == 1 % West wall (x < 0)
            % Filter main dataset
            time_mask = (data(:,1) >= start_date) & (data(:,1) <= end_date);
            region_mask = (x < 0);
            selected_data = data(time_mask & region_mask, :);
            
            % Filter K dataset
            time_mask_K = (data_K(:,1) >= start_date) & (data_K(:,1) <= end_date);
            region_mask_K = (x_K < 0);
            selected_data_K = data_K(time_mask_K & region_mask_K, :);
        else % East wall (x > 0)
            % Filter main dataset
            time_mask = (data(:,1) >= start_date) & (data(:,1) <= end_date);
            region_mask = (x > 0);
            selected_data = data(time_mask & region_mask, :);
            
            % Filter K dataset
            time_mask_K = (data_K(:,1) >= start_date) & (data_K(:,1) <= end_date);
            region_mask_K = (x_K > 0);
            selected_data_K = data_K(time_mask_K & region_mask_K, :);
        end
        
        counts_data = histcounts(selected_data(:,4), bin_edges);
        counts_data_K = histcounts(selected_data_K(:,4), bin_edges);
        
        % Plot both datasets
        b1 = bar(bin_centers, counts_data, 'FaceColor', 'b', 'FaceAlpha', 0.6, 'EdgeColor', 'b');
        hold on;
        b2 = bar(bin_centers, counts_data_K, 'FaceColor', 'r', 'FaceAlpha', 0.6, 'EdgeColor', 'r');
        hold off;
        
        % Set axes with independent limits
        current_max = max([counts_data, counts_data_K]);
        ylim([0 max(current_max*1.1, 1)]); % Ensure at least some height
        xlim([0 2]);
        xticks(0:0.5:2);
        grid on;
        
        % Add title and labels
        title(sprintf('%d', years(y)), 'FontWeight', 'bold');
        if mod(y-1, cols) == 0
            ylabel('Count');
        end
        if y > (rows-1)*cols
            xlabel('Depth (km)');
        end
        
        % Add legend to first subplot only
        if y == 1
            legend([b1, b2], {'William', 'Kaiwen'}, 'Location', 'northwest');
        end
        
        % Add total counts as text
        text(0.1, max(ylim)*0.9, sprintf('W:%d', sum(counts_data)), ...
            'Color', 'b', 'FontSize', 8, 'FontWeight', 'bold');
        text(0.1, max(ylim)*0.8, sprintf('K:%d', sum(counts_data_K)), ...
            'Color', 'r', 'FontSize', 8, 'FontWeight', 'bold');
    end
    
    % Add overall title
    if region == 1
        sgtitle('West Wall (x < 0) - Annual Depth Distribution (2015-2025)', ...
               'FontSize', 14, 'FontWeight', 'bold');
    else
        sgtitle('East Wall (x > 0) - Annual Depth Distribution (2015-2025)', ...
               'FontSize', 14, 'FontWeight', 'bold');
    end
    
    % Adjust subplot spacing
    set(gcf, 'Units', 'normalized');
    h = get(gcf, 'Children');
%    set(h, 'Position', get(h, 'Position').*[1 1 1.05 1.05]);
end