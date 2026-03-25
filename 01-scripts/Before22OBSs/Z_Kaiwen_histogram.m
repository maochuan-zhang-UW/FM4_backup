clc; clear;

% Load the data
%load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/K_W_041525.mat');
load('/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs/A_All/K_W_041525.mat')
data = data_W;

% Define the bin edges for the histogram
dday = 30;
bin_edges = min(data_K(:,1)):dday:max(data_K(:,1));

% Calculate the histogram for Kaiwen
[counts_kaiwen, ~] = histcounts(data_K(:,1), bin_edges);

% Define the three time periods
period1_start = datenum('2015-01-22');
period1_end = datenum('2021-12-31');
period2_start = datenum('2022-01-01');
period2_end = datenum('2022-09-09');
period3_start = datenum('2022-09-10');
period3_end = datenum('2024-09-10');
period4_start = datenum('2024-09-11');

% Plotting
figure; set(gcf, 'Position', [653 595 1158 643]);

% Add background colors for the three periods
axes('Position', [0.13 0.11 0.775 0.815]); % Match the default axes position
% Period 1: 2015-2021 (medium blue)
patch([period1_start period1_end period1_end period1_start], [0 0 20000 20000], ...
      [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on;
% Period 2 and 4: Jan 1, 2022 - Sep 9, 2022 and Sep 11, 2024 onward (medium red)
patch([period2_start period2_end period2_end period2_start], [0 0 20000 20000], ...
      [1 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
patch([period4_start max(data_K(:,1)) max(data_K(:,1)) period4_start], [0 0 20000 20000], ...
      [1 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% Period 3: Sep 10, 2022 - Sep 10, 2024 (medium green)
patch([period3_start period3_end period3_end period3_start], [0 0 20000 20000], ...
      [0.6 1 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Plot Kaiwen's histogram
bar(bin_edges(1:end-1), counts_kaiwen, 'FaceColor', 'r', 'EdgeColor', 'r', ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
ylabel('Event Count');
ylim([0 20000]);
xlim([min(data_K(:,1)) max(data_K(:,1))]);

% Count events in the three periods
events_period1 = sum(data_K(:,1) >= period1_start & data_K(:,1) <= period1_end);
events_period2 = sum(data_K(:,1) >= period2_start & data_K(:,1) <= period2_end) + ...
                sum(data_K(:,1) >= period4_start);
events_period3 = sum(data_K(:,1) >= period3_start & data_K(:,1) <= period3_end);

% Add labels and event counts to the plot
% Period 1 label and count
text(period1_start + (period1_end - period1_start)/2, 18000, ...
     sprintf('Jan 22, 2015 - Dec 31, 2021\nEvents: %d', events_period1), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', ...
     'BackgroundColor', [0.6 0.8 1, 0.5], 'EdgeColor', 'k');

% Period 2 and 4 label and count
text(period2_start + (period3_end - period2_start)/2, 18000, ...
     sprintf('Jan 1, 2022 - April 15, 2025\nEvents: %d', 103595), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', ...
     'BackgroundColor', [1 0.6 0.6, 0.5], 'EdgeColor', 'k');

% Period 3 label and count
text(period3_start + (period3_end - period3_start)/2, 15000, ...
     sprintf('Sep, 2022 - Sep, 2024\nEvents: %d', events_period3), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', ...
     'BackgroundColor', [0.6 1 0.6, 0.5], 'EdgeColor', 'k');

% Add legend below the second period
% legend('Kaiwen', 'Location', 'south', 'Position', ...
%        [0.5, 0.05, 0.1, 0.05], 'HorizontalAlignment', 'center');

% Title and labels
title(['Kaiwen Eqs Count every ' num2str(dday) ' days']);
xlabel('Date');

% Adjust x-axis to display dates
specific_dates = {'2016-01-01', '2017-01-01', '2018-01-01', '2019-01-01', ...
                  '2020-01-01', '2021-01-01', '2022-01-01', '2023-01-01', ...
                  '2024-01-01', '2025-01-01', '2026-01-01'};
date_nums = datenum(specific_dates, 'yyyy-mm-dd');
ax = gca;
ax.XTick = date_nums;
ax.XTickLabel = datestr(ax.XTick, 'yyyy-mm-dd');

% Display the comparison in the console
fprintf('Number of events from Jan 22, 2015, to Dec 31, 2021: %d\n', events_period1);
fprintf('Number of events from Jan 1, 2022 - Sep 9, 2022 and Sep 11, 2024 - onward: %d\n', events_period2);
fprintf('Number of events from Sep 10, 2022, to Sep 10, 2024: %d\n', events_period3);
fprintf('Ratio (Period 2+4 / Period 1): %.2f\n', events_period2 / events_period1);
fprintf('Ratio (Period 3 / Period 1): %.2f\n', events_period3 / events_period1);