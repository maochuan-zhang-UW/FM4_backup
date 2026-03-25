clc; clear;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_all/For_MaleenV3.mat');
data = data_K;

data_K(:,3)=-data_K(:,3);
% Define date limits
startDate = datenum('2022-09-01');
endDate   = datenum('2025-09-01');

% Filter data between 2022 and 2025
mask = data_K(:,1) >= startDate & data_K(:,1) <= endDate;
data_filtered = data_K(mask,:);

% Define the regions
regions.West.Lat = [45.925, 45.955];
regions.West.Lon = [-130.03, -130.006];
regions.East.Lat = [45.93, 45.97];
regions.East.Lon = [-130.00, -129.975];

% Prepare bin edges
dday = 1;
bin_edges = startDate:dday:endDate;

% --- WEST region ---
mask_west = data_filtered(:,2) >= regions.West.Lat(1) & data_filtered(:,2) <= regions.West.Lat(2) & ...
            -data_filtered(:,3) >= regions.West.Lon(1) & -data_filtered(:,3) <= regions.West.Lon(2);
west_data = data_filtered(mask_west,:);
[counts_west, ~] = histcounts(west_data(:,1), bin_edges);

% --- EAST region ---
mask_east = data_filtered(:,2) >= regions.East.Lat(1) & data_filtered(:,2) <= regions.East.Lat(2) & ...
            -data_filtered(:,3) >= regions.East.Lon(1) & -data_filtered(:,3) <= regions.East.Lon(2);
east_data = data_filtered(mask_east,:);
[counts_east, ~] = histcounts(east_data(:,1), bin_edges);

% --- Plotting ---
figure; set(gcf, 'Position', [653 595 1158 643]);

% WEST subplot
subplot(2,1,1);
bar(bin_edges(1:end-1), counts_west, 'FaceColor', 'b', 'EdgeColor', 'b', ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
ylabel('Event Count');
ylim([0 max(counts_west)*1.2]);
xlim([startDate endDate]);
title('West Region Earthquake Count');
xlabel('Date');

% EAST subplot
subplot(2,1,2);
bar(bin_edges(1:end-1), counts_east, 'FaceColor', 'r', 'EdgeColor', 'r', ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
ylabel('Event Count');
ylim([0 max(counts_east)*1.2]);
xlim([startDate endDate]);
title('East Region Earthquake Count');
xlabel('Date');

% Adjust x-axis ticks for both subplots
specific_dates = {'2022-01-01', '2023-01-01', '2024-01-01', '2025-01-01'};
date_nums = datenum(specific_dates, 'yyyy-mm-dd');

ax1 = subplot(2,1,1);
ax1.XTick = date_nums;
ax1.XTickLabel = datestr(ax1.XTick, 'yyyy-mm-dd');
tick_dt = datetime(2022,9,1):calmonths(4):datetime(2025,12,1);
tick_dn = datenum(tick_dt);  % keep using datenum-based axes

% Apply to both subplots
ax1 = subplot(2,1,1);
ax1.XTick = tick_dn;
ax1.XTickLabel = cellstr(datestr(tick_dn,'yyyy-mm'));
ax1.XLim = [startDate endDate];
ax1.XTickLabelRotation = 45;

ax2 = subplot(2,1,2);
ax2.XTick = tick_dn;
ax2.XTickLabel = cellstr(datestr(tick_dn,'yyyy-mm'));
ax2.XLim = [startDate endDate];
ax2.XTickLabelRotation = 45;


% --- Scatter plot of locations for West and East ---
figure; set(gcf, 'Position', [400 400 800 600]);

% West points (blue)
scatter(west_data(:,3), west_data(:,2), 20, 'b', 'filled'); 
hold on;

% East points (red)
scatter(east_data(:,3), east_data(:,2), 20, 'r', 'filled'); 

xlabel('Longitude (°W)');
ylabel('Latitude (°N)');
title('Event Locations: West vs East Regions (2022–2025)');
legend('West Region','East Region','Location','best');
grid on;

% Flip longitude back to negative for display
set(gca, 'XDir','reverse');


% ax2 = subplot(2,1,2);
% ax2.XTick = date_nums;
% ax2.XTickLabel = datestr(ax2.XTick, 'yyyy-mm-dd');
