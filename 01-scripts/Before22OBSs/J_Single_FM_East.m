clear; close all;

% Define paths, dates, and regions
path = '/Users/mczhang/Documents/GitHub/FM4/02-data/';
lonLim = [-130.00 -129.975];
latLim = [45.93 45.97];
regions = struct(...
    'West', struct('Lat', [45.925, 45.955], 'Lon', [-130.03, -130.006]), ...
    'East', struct('Lat', [45.93, 45.97], 'Lon', [-130.00, -129.975]), ...
    'All', struct('Lat', [45.9, 46.01], 'Lon', [-130.00, -129.9]));

% Load caldera rim data
axial_calderaRim;

% Create figure
figure;
set(gcf, 'Position', [280, 378, 600, 500], 'Color', 'white');

% Create axes
ax = axes;
hold on;
box on;
plot(calderaRim(:,1), calderaRim(:,2), '-k', 'linewidth', 3);
Axial_base;

% Load and filter event data
load('/Users/mczhang/Documents/GitHub/FM4/02-data/G_FM/G_All_finalgap.mat');

event = Po_Clu;

% Filter events: mechqual A or B, East region, exclude faultType U, N, R
%event = event(ismember([event.mechqual], {'A', 'B'}));
event([event.lat] < regions.East.Lat(1) | [event.lat] > regions.East.Lat(2)) = [];
event([event.lon] < regions.East.Lon(1) | [event.lon] > regions.East.Lon(2)) = [];
event([event.faultType] == 'U') = [];
event([event.faultType] == 'N') = [];
event([event.faultType] == 'R') = [];


% Plot balloons
radius = 0.0002;
scale = 1.45;
for i = 1:length(event)
    if isempty(event(i).Mw)
        event(i).Mw = 0.1;
    end
    if ~isempty(event(i).avfnorm)
        plot_balloon(event(i).avfnorm, event(i).avslip, event(i).lon, event(i).lat, ...
            event(i).Mw * radius * 1.5, scale, event(i).color2);
    end
end

% Add magnitude legend
i = min(10, length(event)); % Ensure valid index
if i > 0
    plot_balloon(event(i).avfnorm, event(i).avslip, -129.997, 45.935, 2 * radius * 1.5, scale, event(i).color2);
    text(-129.9957, 45.935, 'Mw=2.0');
    plot_balloon(event(i).avfnorm, event(i).avslip, -129.997, 45.934, 1 * radius * 1.5, scale, event(i).color2);
    text(-129.9957, 45.934, 'Mw=1.0');
    plot_balloon(event(i).avfnorm, event(i).avslip, -129.997, 45.933, 0.5 * radius * 1.5, scale, event(i).color2);
    text(-129.9957, 45.933, 'Mw=0.5');
end

% Plot depth points
lon_points = [-130.0000, -129.9998, -130.0000, -129.9998, -129.9998, -129.9977, -129.9934, -129.9894];
lat_points = [45.9521, 45.9546, 45.9575, 45.9614, 45.9660, 45.9682, 45.9687, 45.9689];
depths = [-0.3, -0.6, -0.9, -1.2, -1.5, -1.8, -2.1, -2.4];
for i = 1:length(depths)
    text(lon_points(i), lat_points(i), sprintf('%.1f', depths(i)), ...
        'FontSize', 10, 'Color', 'magenta', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end

% Set axes properties
pbaspect(ax, [diff(lonLim)*cosd(mean(latLim)) diff(latLim) 1]);
axis([-130.001, -129.975, 45.93, 45.97]);
title('East Region All Events 2015 Eruption (Mechqual A or B, Strike-slip)');
xlabel('Longitude', 'fontsize', 14);
ylabel('Latitude', 'fontsize', 14);
set(gca, 'fontSize', 12);
grid on;
colorbar('off');

% Save to PDF
