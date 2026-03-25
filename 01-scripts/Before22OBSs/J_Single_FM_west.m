clear; close all;

% Define paths, dates, and regions
path = '/Users/mczhang/Documents/GitHub/FM4/02-data/';
lonLim = [-130.029, -130.008];
latLim = [45.93, 45.953];
regions = struct('West', struct('Lat', latLim, 'Lon', lonLim));

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

% Load and filter event data
load('/Users/mczhang/Documents/GitHub/FM4/02-data/G_FM/G_All_finalgap.mat');
event = Po_Clu;

% Filter events: mechqual A or B, within West region, exclude faultType U
event = event(ismember([event.mechqual], {'A', 'B'}));
event([event.lat] < latLim(1) | [event.lat] > latLim(2)) = [];
event([event.lon] < lonLim(1) | [event.lon] > lonLim(2)) = [];
event([event.faultType] == 'U') = [];

% Sort events by depth
[~, I] = sort([event.depth], 'ascend');
event = event(I);

% Common settings for region limits
x1 = lonLim(1);
x2 = lonLim(2);
y1 = latLim(1);
y2 = latLim(2);
radius = 0.0002;

% Plot balloons
for i = 1:length(event)
    if isempty(event(i).Mw)
        event(i).Mw = 0.1;
    end
    plot_balloon(event(i).avfnorm, event(i).avslip, event(i).lon, event(i).lat, ...
        event(i).Mw * radius * 1.5, 1.0, event(i).color);
end

% Add magnitude legend
i = 10; % Use event 10 for consistent color
plot_balloon(event(i).avfnorm, event(i).avslip, -130.012, 45.937, 2 * radius * 1.5, 1.0, event(i).color);
text(-130.0115, 45.937, 'Mw=2.0');
plot_balloon(event(i).avfnorm, event(i).avslip, -130.012, 45.936, 1 * radius * 1.5, 1.0, event(i).color);
text(-130.0115, 45.936, 'Mw=1.0');
plot_balloon(event(i).avfnorm, event(i).avslip, -130.012, 45.935, 0.5 * radius * 1.5, 1.0, event(i).color);
text(-130.0115, 45.935, 'Mw=0.5');

% Plot fissure line
x = [-130.0308; -130.0056];
y = [45.944; 45.944];
plot(x, y, 'r-.', 'LineWidth', 2);

% Set axes properties
axis equal;
xlim([x1, x2]);
ylim([y1, y2]);
title('West Region All Events 2015 Eruption (Mechqual A or B)');
xlabel('Longitude', 'fontsize', 14);
ylabel('Latitude', 'fontsize', 14);
specificYTicks = [45.930, 45.935, 45.940, 45.945, 45.950];
yticks(specificYTicks);
specificTicks = [-130.025, -130.02, -130.015, -130.01];
xticks(specificTicks);
set(gca, 'fontSize', 11);
grid on;
colorbar('off');