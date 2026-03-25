clc; clear; close all;

% Load W1 and W2 clusters
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_ID/A_W1.mat');
W1_events = Felix;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_ID/A_W2.mat');
W2_events = Felix;

% Combine all W1 and W2 events
all_events = [W1_events, W2_events];

% Convert to Cartesian coordinates
[x_all, y_all] = latlon2xy([all_events.lat], [all_events.lon]);
XY_all = [x_all', y_all'];

% Define the boundary line (your provided coordinates)
boundary_start = [-0.126, -1.922];  % x1,y1
boundary_end = [-6.376, 1.323];     % x2,y2

% Calculate distance of each point to the boundary line
boundary_vector = boundary_end - boundary_start;
boundary_length = norm(boundary_vector);
unit_vector = boundary_vector / boundary_length;

% Vector from start point to each event
vectors_to_events = XY_all - boundary_start;

% Projection lengths (scalar projections)
proj_lengths = vectors_to_events * unit_vector';

% Clamp projections to line segment boundaries
proj_lengths = max(0, min(boundary_length, proj_lengths));

% Find closest points on line segment
closest_points = boundary_start + proj_lengths * unit_vector;

% Calculate perpendicular distances
distances = sqrt(sum((XY_all - closest_points).^2, 2));

% Find events within 0.5 km of boundary
is_near_boundary = distances < 0.5;
boundary_events = all_events(is_near_boundary);
clear Felix;
Felix=boundary_events;
% Save results
save('/Users/mczhang/Documents/GitHub/FM4/02-data/A_ID/A_W12.mat', 'Felix');
fprintf('Found %d events within 0.5 km of W1-W2 boundary\n', length(boundary_events));

% Visualization
figure;
set(gcf, 'position', [100, 100, 800, 800]);
hold on;

% Plot all events
scatter(XY_all(:,1), XY_all(:,2), 20, 'k', 'filled', 'DisplayName', 'All Events');

% Plot boundary line
plot([boundary_start(1), boundary_end(1)], [boundary_start(2), boundary_end(2)], ...
    'r-', 'LineWidth', 2, 'DisplayName', 'W1-W2 Boundary');

% Plot boundary events
scatter(XY_all(is_near_boundary,1), XY_all(is_near_boundary,2), 40, 'g', 'd', ...
    'filled', 'DisplayName', 'Boundary Events');

% Add reference elements
axial_calderaRim;
[calderaRim(:,2), calderaRim(:,1)] = latlon2xy(calderaRim(:,2), calderaRim(:,1));
plot(calderaRim(:,2), calderaRim(:,1), 'k', 'LineWidth', 2);

sta = axial_stationsNewOrder;
sta = sta(1:7); 
for i = 1:length(sta)
    [sta(i).x, sta(i).y] = latlon2xy_no_rotate([sta(i).lat], [sta(i).lon]);
    plot(sta(i).x, sta(i).y, 's', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
    text(sta(i).x, sta(i).y, sta(i).name(3:end), 'Color', 'k');
end

title('W1-W2 Boundary Events Identification');
legend('Location', 'best');
axis equal;
xlim([min(XY_all(:,1))-0.5, max(XY_all(:,1))+0.5]);
ylim([min(XY_all(:,2))-0.5, max(XY_all(:,2))+0.5]);
grid on;
set(gca, 'FontSize', 12);
xlabel('X km');
ylabel('Y km');
hold off;