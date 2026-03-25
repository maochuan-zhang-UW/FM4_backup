clc; clear; close all;


lonLim = [-130.031 -129.97];
latLim = [45.92 45.972];

%% ===========================================
% PART 1 — LOAD 1D/3D KAGAN DATA
%% ===========================================
load('/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs/G_FM/G_HASH_All.mat');

event = event1;
event([event.mechqual] == 'D' | [event.mechqual] == 'C') = [];
clear event1;

load('/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs/G_FM/G_HASH_All_ML_sameClusterasbefore.mat');

%% ---- Compute REAL Kagan angles (1-D vs 3-D) ----
for i = 1:length(event)
    ind = find([event1.id] == event(i).id );
    if ~isempty(ind) && ind(1) <= 4090
        event(i).kg = kagan( ...
            [event(i).avmech(1), event(i).avmech(2), event(i).avmech(3)], ...
            [event1(ind(1)).avmech(1), event1(ind(1)).avmech(2), event1(ind(1)).avmech(3)]);
    end
end

% Remove empty kg
event(arrayfun(@(x) isempty(x.kg), event)) = [];

%% ---- Extract data ----
kg  = real([event.kg]);
lat = [event.lat];
lon = [event.lon];

%% ===========================================
% PART 2 — HISTOGRAM WITH MEAN & MEDIAN
%% ===========================================
kg_mean   = mean(kg);
kg_median = median(kg);

figure('color','w','position',[200 200 1200 500]);

subplot(1,2,1)
histogram(kg,20,'FaceColor',[0.2 0.4 0.8],'EdgeColor','k');
xlabel('Kagan Angle (°)');
ylabel('Count');
title(sprintf('Kagan Angle Distribution  |  Mean = %.2f°   Median = %.2f°', ...
    kg_mean, kg_median), 'FontWeight','bold');
grid on; box on;

%% ===========================================
% PART 3 — MAP WITH CALDERA + KAGAN COLOR
%% ===========================================
subplot(1,2,2)
hold on

% ---- Plot caldera rim FIRST ----
axial_calderaRim;
plot(calderaRim(:,1), calderaRim(:,2), '-k', 'linewidth', 3, 'HandleVisibility', 'off');

% ---- Plot earthquake locations colored by Kagan ----
scatter(lon, lat, 25, kg, 'filled');

colormap(jet);
colorbar;
caxis([0 90]);

xlabel('Longitude');
ylabel('Latitude');
title('Event Location Colored by Kagan Angle');

axis([-130.0300 -129.97 45.9200 45.9700]);
axis equal;
grid on; box on;
set(gca,'FontSize',12,'LineWidth',1.2)