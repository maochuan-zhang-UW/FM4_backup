%% Main Script: One Figure with All Events
clear; close all;

% Define cutoff dates and regions
date_BF = datenum(2025, 4, 24, 8, 0, 0);
date_DR = datenum(2015, 5, 19);
regions = struct(...
    'West', struct('Lat', [45.930, 45.953], 'Lon', [-130.029, -130.008]), ...
    'East', struct('Lat', [45.971, 45.930], 'Lon', [-130.0015, -129.975]), ...
    'ID', struct('Lat', [45.921, 45.929], 'Lon', [-130.004, -129.975]));
% lonLim = [-130.04 -129.97];
% latLim = [45.908 46.001];
% Define latitude and longitude for each station
latitudes = [45.9336, 45.9338, 45.9547, 45.9496, 45.9397, 45.9361, 45.9257, ...
             46.01933, 46.000625, 45.992922, 45.988743, 45.97439, 45.982938, 45.969067, ...
             45.964147, 45.97084, 45.960338, 45.94959, 45.915367, 45.907817, 45.899965];

longitudes = [-129.9992, -130.0141, -130.0089, -129.9797, -129.9738, -129.9785, -129.9780, ...
              -130.005375, -130.022257, -129.992782, -130.048395, -129.977972, -130.014128, -130.029852, ...
              -130.004503, -130.062393, -129.949987, -130.033977, -130.021773, -129.971555, -130.007573];


lonLim = [min(longitudes) max(longitudes)];
latLim = [min(latitudes) max(latitudes)];

% Load the event data
load('/Users/mczhang/Documents/GitHub/FM4/02-data/G_FM/G_3D.mat');
%load('/Users/mczhang/Documents/GitHub/FM4/02-data/G_FM/G_FM_SKHASHV2.mat')
%event1D = event1;
%load('/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs/G_FM/G_HASH_All.mat');
event1D = event1;

% load('/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs/G_FM/G_All_finalgap.mat');
% event1D = Po_Clu;
date1=datenum(2022, 9, 1, 0, 0, 0);
date2=datenum(2023, 9, 1, 0, 0, 0);
%event1D([event1D.on] <= date1)=[];
%event1D([event1D.on] >= date2)=[];
event1D([event1D.faultType] =='U')=[];

% Filter events by mechquality 'A' or 'B' and assign colors'B','C','D'
event1D = event1D(ismember({event1D.mechqual}, {'A','B'}));
event1D = event1D(ismember({event1D.faultType}, {'N','R','S'}));
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
set(ax, 'XTick', -130.06:0.01:-129.94);
set(ax, 'GridLineStyle', '-', 'LineWidth', 0.5, 'GridColor', [0.5 0.5 0.5]);
grid(ax, 'on');

% Filter events spatially
events = event1D([event1D.lat] <= 46 & [event1D.lon] >= -130.04);
nmec = length(events);

% Plot each event
radius = 0.0005;
scale_event = 1.3;
for i = 1:nmec
    if ~isempty(events(i).avfnorm)
        plot_balloon(events(i).avfnorm, events(i).avslip, ...
            events(i).lon, events(i).lat, radius*2, scale_event, events(i).color3);
        hold(ax, 'on');
    end
    %text(events(i).lon+0.001,events(i).lat,num2str(events(i).id));
    %text(events(i).lon+0.006,events(i).lat,num2str(events(i).depth,2));

end

hold(ax, 'on');
% Define station codes
station_codes = {'AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1', ...
                 '01A', '02A', '03A', '04A', '05A', '06A', '07A', ...
                 '08A', '09A', '10A', '11A', '12A', '13A', '14A'};

% Plot stations: first 7 in full black, rest in lighter gray
plot(longitudes(1:7), latitudes(1:7), 'sk', ...
    'MarkerFaceColor', 'black', 'MarkerSize', 14);
plot(longitudes(8:end), latitudes(8:end), 'sk', ...
    'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerSize', 14);

%plot(longitudes, latitudes, 'sk', ...
%  'markerfacecolor', 'black', 'markersize', 14);
% Set labels and title
xlabel(ax, 'Longitude', 'fontsize', 14);
ylabel(ax, 'Latitude', 'fontsize', 14);
title(ax, [num2str(nmec) ' FMs (2022.9-2023.9, Mechqual A or B)'], 'fontsize', 14);

% % Add legend
% legendX = lonLim(1) + 0.008;
% legendYs = linspace(latLim(2)-0.002, latLim(2)-0.01, 4);
% faultTypes = {'N', 'R', 'S', 'U'};
% faultLabels = {'N - Normal', 'R - Reverse', 'S - Strike-slip', 'U - Undefined'};
% for k = 1:length(faultTypes)
%     idx = find(strcmp({events.faultType}, faultTypes{k}), 1, 'first');
%     if ~isempty(idx)
%         repEvent = events(idx);
%         plot_balloon(repEvent.avfnorm, repEvent.avslip, ...
%             legendX, legendYs(k), radius, scale_event, repEvent.color2);
%     else
%         defaultColors = containers.Map(...
%             {'N', 'R', 'S', 'U'}, {[1,0,0], [0,0,1], [0,1,0], [0,0,0]});
%         plot_balloon(0, 0, legendX, legendYs(k), radius, scale_event, defaultColors(faultTypes{k}));
%     end
%     text(legendX + 0.003, legendYs(k), faultLabels{k}, ...
%         'HorizontalAlignment', 'left', 'FontSize', 10);
% end

% % Draw legend rectangle
% legendLeft = legendX - 0.002;
% legendRight = legendX + 0.020;
% legendBottom = min(legendYs) - 0.002;
% legendTop = max(legendYs) + 0.001;
% legendWidth = legendRight - legendLeft;
% legendHeight = legendTop - legendBottom;
% rectangle(ax, 'Position', [legendLeft, legendBottom, legendWidth, legendHeight], ...
%     'EdgeColor', 'k', 'LineWidth', 0.5);

% Plot regions
% regionNames = fieldnames(regions);
% for i = 1:length(regionNames)
%     region = regions.(regionNames{i});
%     latitudes = [region.Lat(1), region.Lat(1), region.Lat(2), region.Lat(2), region.Lat(1)];
%     longitudes = [region.Lon(1), region.Lon(2), region.Lon(2), region.Lon(1), region.Lon(1)];
%     plot3(longitudes, latitudes, zeros(size(latitudes)), '--', 'LineWidth', 2, ...
%         'DisplayName', [regionNames{i} ' Region']);
% end