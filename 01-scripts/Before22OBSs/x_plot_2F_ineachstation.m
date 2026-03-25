clear; clc;close all

% Load the saved file
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/A_Kaiwen_phase2F.mat', 'Felix');

% Define the station names and their coordinates
stations = {
    struct('name', 'AXCC1', 'lat', 45.9547, 'lon', -130.0089, 'code', 'CC1'),
    struct('name', 'AXEC1', 'lat', 45.9496, 'lon', -129.9797, 'code', 'EC1'),
    struct('name', 'AXEC2', 'lat', 45.9397, 'lon', -129.9738, 'code', 'EC2'),
    struct('name', 'AXEC3', 'lat', 45.9361, 'lon', -129.9785, 'code', 'EC3'),
    struct('name', 'AXAS1', 'lat', 45.9336, 'lon', -129.9992, 'code', 'AS1'),
    struct('name', 'AXAS2', 'lat', 45.9338, 'lon', -130.0141, 'code', 'AS2'),
    struct('name', 'AXID1', 'lat', 45.9257, 'lon', -129.9780, 'code', 'ID1'),
    struct('name', 'AX01A', 'lat', 46.01933, 'lon', -130.005375, 'code', '01A'),
    struct('name', 'AX02A', 'lat', 46.000625, 'lon', -130.029257, 'code', '02A'),
    struct('name', 'AX03A', 'lat', 45.992922, 'lon', -129.992782, 'code', '03A'),
    struct('name', 'AX04A', 'lat', 45.988743, 'lon', -130.048395, 'code', '04A'),
    struct('name', 'AX05A', 'lat', 45.97439, 'lon', -129.977972, 'code', '05A'),
    struct('name', 'AX06A', 'lat', 45.982938, 'lon', -130.014128, 'code', '06A'),
    struct('name', 'AX07A', 'lat', 45.969067, 'lon', -130.029852, 'code', '07A'),
    struct('name', 'AX08A', 'lat', 45.964147, 'lon', -130.004503, 'code', '08A'),
    struct('name', 'AX09A', 'lat', 45.97084, 'lon', -130.062393, 'code', '09A'),
    struct('name', 'AX10A', 'lat', 45.960338, 'lon', -129.949987, 'code', '10A'),
    struct('name', 'AX11A', 'lat', 45.94959, 'lon', -130.039377, 'code', '11A'),
    struct('name', 'AX12A', 'lat', 45.915367, 'lon', -130.021777, 'code', '12A'),
    struct('name', 'AX13A', 'lat', 45.907817, 'lon', -129.971555, 'code', '13A'),
    struct('name', 'AX14A', 'lat', 45.899965, 'lon', -130.007573, 'code', '14A'),
    struct('name', 'AX15A', 'lat', 45.91972, 'lon', -129.93942, 'code', '15A')
};

% Create a figure
figure;
set(gcf, 'Position', [100, 100, 2000, 1500]); % Set figure size

% Loop through each station and create a subplot
for i = 1:length(stations)
    subplot(5, 5, i); % Create subplot in a 5x5 grid
    hold on;
    
    % Extract the station information
    station = stations{i};
    stationCode = station.code;
    
    % Initialize empty arrays for latitude and longitude
    lat = [];
    lon = [];
    
    % Loop through events and collect those with non-zero DDt
    for j = 1:length(Felix)
        if isfield(Felix(j), ['DDt_' stationCode]) && Felix(j).(['DDt_' stationCode]) ~= 0
            lat = [lat; Felix(j).lat];
            lon = [lon; Felix(j).lon];
        end
    end
    
    % Scatter plot the events
    scatter(lon, lat, 5, 'filled');
    plot(station.lon, station.lat, 'r*', 'MarkerSize', 6); % Plot the station with a red star
    
    % Add grid, labels, and title
    grid on;
    xlabel('Longitude');
    ylabel('Latitude');
    title(sprintf('%s (N=%d)', station.name, length(lat)));
    hold on;
    axial_calderaRim;
    plot(calderaRim(:,1), calderaRim(:,2), 'k', 'LineWidth', 3);
end

% Adjust layout for better visualization
sgtitle('Events with Non-zero DDt for Each Station');
