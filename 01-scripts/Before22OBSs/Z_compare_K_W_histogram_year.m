clear; close all;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/K_W_041525.mat');
data = data_W;

% Define parameters
bin_edges = 0:0.1:2;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;
years = 2015:2025;
num_years = length(years);

% Create figure with subplots
figure;
set(gcf, 'Position', [100, 100, 1200, 900]);

% Determine subplot arrangement (3x4 grid for 11 years)
rows = 3;
cols = 4;

% Plot each year's data with independent y-axis
for y = 1:num_years
    subplot(rows, cols, y);
    
    % Get data for current year
    start_date = datenum(years(y), 1, 1);
    end_date = datenum(years(y), 12, 31);
    
    selected_data = data(data(:,1) >= start_date & data(:,1) <= end_date, :);
    selected_data_K = data_K(data_K(:,1) >= start_date & data_K(:,1) <= end_date, :);
    
    counts_data = histcounts(selected_data(:,4), bin_edges);
    counts_data_K = histcounts(selected_data_K(:,4), bin_edges);
    
    % Plot both datasets
    b1 = bar(bin_centers, counts_data, 'FaceColor', 'b', 'FaceAlpha', 0.6, 'EdgeColor', 'b');
    hold on;
    b2 = bar(bin_centers, counts_data_K, 'FaceColor', 'r', 'FaceAlpha', 0.6, 'EdgeColor', 'r');
    hold off;
    
    % Set axes - let each subplot determine its own y-limits
    current_max = max([counts_data, counts_data_K]);
    ylim([0 current_max*1.1]);
    xlim([0 2]);
    xticks(0:0.5:2);
    grid on;
    
    % Add title and labels
    title(sprintf('%d', years(y)), 'FontWeight', 'bold');
    if mod(y-1, cols) == 0  % Only show ylabel for leftmost plots
        ylabel('Count');
    end
    if y > (rows-1)*cols  % Only show xlabel for bottom plots
        xlabel('Depth (km)');
    end
    
    % Add legend to first subplot only
    if y == 1
        legend([b1, b2], {'William', 'Kaiwen'}, 'Location', 'northwest');
    end
    
    % Add total counts as text
    text(0.1, current_max*0.9, sprintf('W:%d', sum(counts_data)), ...
        'Color', 'b', 'FontSize', 8, 'FontWeight', 'bold');
    text(0.1, current_max*0.8, sprintf('K:%d', sum(counts_data_K)), ...
        'Color', 'r', 'FontSize', 8, 'FontWeight', 'bold');
end

% Add overall title
sgtitle('Annual Depth Distribution Comparison (2015-2025) - Independent Y-Axis', ...
       'FontSize', 14, 'FontWeight', 'bold');

% Adjust subplot spacing
set(gcf, 'Units', 'normalized');
h = get(gcf, 'Children');
%set(h, 'Position', get(h, 'Position').*[1 1 1.05 1.05]);  % Add small margin