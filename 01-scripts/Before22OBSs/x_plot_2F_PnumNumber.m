clear; clc;

% Load the saved file
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/A_Kaiwen_phase2F.mat', 'Felix');

% Define the range of Pnum thresholds (5 to 21)
thresholds = 5:21;

% Create a figure
figure;
set(gcf, 'Position', [100, 100, 1800, 1400]); % Set figure size

% Loop through each threshold and create a subplot
for i = 1:length(thresholds)
    subplot(4, 5, i); % Create subplot in a 4x5 grid
    hold on;
    
    % Initialize empty arrays for latitude and longitude
    lat = [];
    lon = [];
    dep = [];
    
    % Loop through events and collect those with Pnum > threshold
    for j = 1:length(Felix)
        if isfield(Felix(j), 'Pnum') && Felix(j).Pnum > thresholds(i)
            lat = [lat; Felix(j).lat];
            lon = [lon; Felix(j).lon];
            dep = [dep; Felix(j).depth];
        end
    end
    
    % Scatter plot the events
    scatter3(lon, lat,dep, 5, 'filled');
    
    hold on;
    axial_calderaRim;
    plot(calderaRim(:,1), calderaRim(:,2), 'k', 'LineWidth', 3);
    grid on;
    view(0,90);
    xlabel('Longitude');
    ylabel('Latitude');
    title(sprintf('Pnum > %d (N=%d)', thresholds(i), length(lat)));
end

% Adjust layout for better visualization
clear
sgtitle('Events with Pnum Greater Than Threshold');
