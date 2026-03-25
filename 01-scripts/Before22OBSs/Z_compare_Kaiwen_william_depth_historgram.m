clear; close all;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/K_W_041525.mat');
data = data_W;

% Set your date range (modify these as needed)
start_date = datenum('2022-01-01');
end_date = datenum('2024-01-01');

% Filter both datasets for the time period
selected_data = data(data(:,1) >= start_date & data(:,1) <= end_date, :);
selected_data_K = data_K(data_K(:,1) >= start_date & data_K(:,1) <= end_date, :);

% Define bin edges for depth (0:0.1:2)
bin_edges = 0:0.1:2;
bin_centers = bin_edges(1:end-1) + diff(bin_edges)/2;

% Calculate histograms
counts_data = histcounts(selected_data(:,4), bin_edges);
counts_data_K = histcounts(selected_data_K(:,4), bin_edges);

% Create figure
figure;
set(gcf, 'Position', [100, 100, 900, 600]);

% Find maximum value for consistent y-axis scaling
max_count = max([counts_data, counts_data_K]);

% Plot both datasets on the same y-axis
bar(bin_centers, counts_data, 'FaceColor', 'b', 'FaceAlpha', 0.6, 'EdgeColor', 'b', 'BarWidth', 1);
hold on;
bar(bin_centers, counts_data_K, 'FaceColor', 'r', 'FaceAlpha', 0.6, 'EdgeColor', 'r', 'BarWidth', 1);
hold off;

% Create dynamic title with formatted dates
start_str = datestr(start_date, 'yyyy.mm.dd');
end_str = datestr(end_date, 'yyyy.mm.dd');
title_str = sprintf('Depth Distribution Comparison (%s - %s)', start_str, end_str);

% Plot settings
ylabel('Number of Events', 'FontWeight', 'bold');
ylim([0 max_count*1.1]);
xlabel('Depth (km)', 'FontWeight', 'bold');
title(title_str, 'FontWeight', 'bold');
grid on;
xlim([0 2]);
xticks(0:0.2:2);

% Add legend
legend({'William Data', 'Kaiwen Data'}, 'Location', 'northeast');

% Display total counts
total_william = sum(counts_data);
total_kaiwen = sum(counts_data_K);
text(0.1, max_count*0.95, sprintf('William total: %d', total_william), ...
    'Color', 'b', 'FontSize', 10, 'FontWeight', 'bold');
text(0.1, max_count*0.85, sprintf('Kaiwen total: %d', total_kaiwen), ...
    'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');