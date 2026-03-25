%% Main Script: One Figure with All Events
clear; close all;
% Define cutoff dates and regions
date_BF = datenum(2025, 4, 24, 8, 0, 0);
date_DR = datenum(2015, 5, 19);
regions = struct(...
    'West', struct('Lat', [45.930, 45.953], 'Lon', [-130.029, -130.008]), ...
    'East', struct('Lat', [45.971, 45.930], 'Lon', [-130.0015, -129.975]), ...
    'ID', struct('Lat', [45.921, 45.929], 'Lon', [-130.004, -129.975]));
lonLim = [-130.031 -129.97];
latLim = [45.92 45.972];

% Load the event data
load('/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs/G_FM/G_HASH_All.mat');
%load('/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs/G_FM/G_HASH_All_ML_sameClusterasbefore.mat');

event1D = event1;
% event1D = Po_Clu;
% date1=datenum(2022, 9, 1, 0, 0, 0);
% date2=datenum(2023, 9, 1, 0, 0, 0);
% event1D([event1D.on] <= date1)=[];
% event1D([event1D.on] >= date2)=[];
event1D([event1D.faultType] =='U')=[];

% Filter events by mechquality 'A' or 'B' and assign colors
event1D = event1D(ismember({event1D.mechqual}, {'A', 'B'}));
for i = 1:length(event1D)
    event1D(i).color3 = event1D(i).color2;
end

% Create figure
figure('Position', [744, 309, 786, 914], ...
    'InvertHardcopy', 'off', ...
    'Color', 'white');

% Create axes
ax = axes;
basemap_2015v2(lonLim, latLim, 100, [0 0], 1, false, ax);
pbaspect(ax, [diff(lonLim)*cosd(mean(latLim)) diff(latLim) 1]);
hold(ax, 'on');
set(ax, 'XTick', -130.03:0.01:-129.97);
set(ax, 'GridLineStyle', '-', 'LineWidth', 0.5, 'GridColor', [0.5 0.5 0.5]);
grid(ax, 'on');

% Filter events spatially
events = event1D([event1D.lat] <= 45.97 & [event1D.lon] >= -130.03);
nmec = length(events);

% Plot each event
radius = 0.0005;
scale_event = 1.3;
for i = 1:nmec
    if ~isempty(events(i).avfnorm)
        plot_balloon(events(i).avfnorm, events(i).avslip, ...
            events(i).lon, events(i).lat, radius, scale_event, events(i).color3);
        hold(ax, 'on');
    end
end

% Set labels and title
xlabel(ax, 'Longitude', 'fontsize', 14);
ylabel(ax, 'Latitude', 'fontsize', 14);
title(ax, [num2str(nmec) ' FMs (2022.1-2025.3, Mechqual A or B)'], 'fontsize', 14);

% Add legend
legendX = lonLim(1) + 0.008;
legendYs = linspace(latLim(2)-0.002, latLim(2)-0.01, 4);
faultTypes = {'N', 'R', 'S', 'U'};
faultLabels = {'N - Normal', 'R - Reverse', 'S - Strike-slip', 'U - Undefined'};
for k = 1:length(faultTypes)-1
    idx = find(strcmp({events.faultType}, faultTypes{k}), 1, 'first');
    if ~isempty(idx)
        repEvent = events(idx);
        plot_balloon(repEvent.avfnorm, repEvent.avslip, ...
            legendX, legendYs(k), radius*2, scale_event, repEvent.color2);
    else
        defaultColors = containers.Map(...
            {'N', 'R', 'S', 'U'}, {[1,0,0], [0,0,1], [0,1,0], [0,0,0]});
        plot_balloon(0, 0, legendX, legendYs(k), radius*2, scale_event, defaultColors(faultTypes{k}));
    end
    text(legendX + 0.003, legendYs(k), faultLabels{k}, ...
        'HorizontalAlignment', 'left', 'FontSize', 22);
end

% Draw legend rectangle
legendLeft = legendX - 0.002;
legendRight = legendX + 0.018;
legendBottom = min(legendYs) ;%- 0.002;
legendTop = max(legendYs) + 0.0015;
legendWidth = legendRight - legendLeft;
legendHeight = legendTop - legendBottom;
rectangle(ax, 'Position', [legendLeft, legendBottom, legendWidth, legendHeight], ...
    'EdgeColor', 'k', 'LineWidth', 0.5);

% Plot regions
regionNames = fieldnames(regions);
for i = 1:length(regionNames)
    region = regions.(regionNames{i});
    latitudes = [region.Lat(1), region.Lat(1), region.Lat(2), region.Lat(2), region.Lat(1)];
    longitudes = [region.Lon(1), region.Lon(2), region.Lon(2), region.Lon(1), region.Lon(1)];
    plot3(longitudes, latitudes, zeros(size(latitudes)), '--', 'LineWidth', 2, ...
        'DisplayName', [regionNames{i} ' Region']);
end
set(gca, "FontSize",22);