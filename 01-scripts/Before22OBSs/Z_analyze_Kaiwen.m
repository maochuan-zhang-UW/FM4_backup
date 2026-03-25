% Load or assume Felix struct is already in workspace
% Felix is a 1x90397 struct array
clear;close all;clc;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_Kaiwen.mat');
%% 1. Time Distribution Histogram
% Extract timestamps
time_data = [Felix.on];

% Check if time_data is in datenum format; convert to datetime if needed
if max(time_data) > 1e5  % Rough check for datenum
    time_datetime = datetime(time_data, 'ConvertFrom', 'datenum');
else
    time_datetime = datetime(time_data); % Adjust based on actual format
end

% Plot histogram of time distribution
figure('Name', 'Time Distribution');
histogram(time_datetime, 'BinMethod', 'auto');
xlabel('Time');
ylabel('Number of Events');
title('Histogram of Event Times');
grid on;

%% 2. Spatial Distribution (Lat, Lon, Depth)
lat = [Felix.lat];
lon = [Felix.lon];
depth = [Felix.depth];

% Scatter plot for lat-lon
figure('Name', 'Spatial Distribution');
subplot(2,1,1);
scatter(lon, lat, 5, 'filled');
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Geographical Distribution of Events');
grid on;

% Histogram for depth
subplot(2,1,2);
histogram(depth, 'BinMethod', 'auto');
xlabel('Depth (km)');
ylabel('Number of Events');
title('Depth Distribution');
grid on;

%% 3. Pnum and Snum Analysis
Pnum = [Felix.Pnum];
Snum = [Felix.Snum];

% Find events where Pnum or Snum > 5
idx_Pnum = Pnum > 4;
idx_Snum = Snum > 4;
idx_either = idx_Pnum & idx_Snum;

% Summary statistics
num_Pnum_gt5 = sum(idx_Pnum);
num_Snum_gt5 = sum(idx_Snum);
num_either_gt5 = sum(idx_either);

fprintf('Number of events with Pnum > 5: %d\n', num_Pnum_gt5);
fprintf('Number of events with Snum > 5: %d\n', num_Snum_gt5);
fprintf('Number of events with either Pnum > 5 or Snum > 5: %d\n', num_either_gt5);

% Optional: List some example events
if num_either_gt5 > 0
    fprintf('\nExample events with Pnum > 5 or Snum > 5:\n');
    fprintf('ID\tTime\t\t\tLat\tLon\tDepth\tPnum\tSnum\n');
    example_idx = find(idx_either, 5, 'first'); % Show up to 5 examples
    for i = example_idx
        fprintf('%s\t%s\t%.2f\t%.2f\t%.2f\t%d\t%d\n', ...
            Felix(i).ID, datestr(Felix(i).on), ...
            Felix(i).lat, Felix(i).lon, Felix(i).depth, ...
            Felix(i).Pnum, Felix(i).Snum);
    end
end

%% 4. Optional: Save figures
% Assume Felix struct is in workspace
% Felix is a 1x90397 struct array with fields lat, lon, depth

%% Extract location data
lat = [Felix.lat];
lon = [Felix.lon];
depth = [Felix.depth];

%% Create 3D Scatter Plot
figure('Name', '3D Location of Events');
scatter3(lon, lat,-depth, 1, 'filled', 'MarkerFaceColor', 'b');
xlabel('Longitude (°)');
ylabel('Latitude (°)');
zlabel('Depth (km)');
title('3D Spatial Distribution of Events');
grid on;
%axis equal; % Equal scaling for better visualization
view(-45, 30); % Adjust view angle for better perspective

% Optional: Invert z-axis if depth is positive downward
set(gca, 'ZDir', 'reverse'); % Depth increases downward

%% Optional: Save figure
% Uncomment to save
% saveas(gcf, '3D_Location_Plot.png');
% Uncomment to save
% saveas(figure(1), 'Time_Distribution.png');
% saveas(figure(2), 'Spatial_Distribution.png');