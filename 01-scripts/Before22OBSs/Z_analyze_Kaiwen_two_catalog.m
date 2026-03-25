clc; clear; close all;

% Load the data
load('/Users/mczhang/Documents/GitHub/FM4/02-data/Kaiwen_eventData_website.mat')
Felix_DD = eventData; clear Felix;

load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_wavelarge5.mat');
Felix_Wa = Felix; clear Felix;

% Find common IDs
commonIDs = intersect([Felix_DD.eventID], [Felix_Wa.ID]);

if isempty(commonIDs)
    error('No common IDs found between Felix_DD and Felix_Wa');
end

% Get indices of common events in both datasets
[~, idx_DD] = ismember(commonIDs, [Felix_DD.eventID]);
[~, idx_Wa] = ismember(commonIDs, [Felix_Wa.ID]);

% Create a copy of Felix_Wa before modification
Felix_Wa_original = Felix_Wa;

% STEP 1: Replace coordinates in Felix_Wa with values from Felix_DD for common events
for i = 1:length(idx_Wa)
    Felix_Wa(idx_Wa(i)).lat = Felix_DD(idx_DD(i)).lat;
    Felix_Wa(idx_Wa(i)).lon = Felix_DD(idx_DD(i)).lon;
    Felix_Wa(idx_Wa(i)).depth = Felix_DD(idx_DD(i)).depth;
end

% STEP 2: Remove events from Felix_Wa that don't have common IDs
% Create logical index of events to keep (those with common IDs)
keep_mask = ismember([Felix_Wa.ID], commonIDs);

% Apply the filter
Felix_Wa_filtered = Felix_Wa(keep_mask);

% Verification
fprintf('Original Felix_Wa had %d events\n', length(Felix_Wa));
fprintf('Filtered Felix_Wa has %d events (only common IDs kept)\n', length(Felix_Wa_filtered));
fprintf('%d events were removed\n', length(Felix_Wa) - length(Felix_Wa_filtered));

% Plot comparison
figure;

% Before processing
subplot(1, 3, 1);
[x_Wa_orig, y_Wa_orig] = latlon2xy([Felix_Wa_original.lat], [Felix_Wa_original.lon]);
[x_DD_plot, y_DD_plot] = latlon2xy([Felix_DD.lat], [Felix_DD.lon]);

scatter(x_Wa_orig, y_Wa_orig, 5, 'b'); hold on;
scatter(x_DD_plot, y_DD_plot, 5, 'r');
scatter(x_Wa_orig(idx_Wa), y_Wa_orig(idx_Wa), 50, 'k', 'filled');
title('Original Data');
legend('Felix\_Wa', 'Felix\_DD', 'Common IDs');
axis equal; grid on;
xlim([-3 3]); ylim([-4.5 4.5]);

% After coordinate replacement
subplot(1, 3, 2);
[x_Wa_new, y_Wa_new] = latlon2xy([Felix_Wa.lat], [Felix_Wa.lon]);

scatter(x_Wa_new, y_Wa_new, 5, 'b'); hold on;
scatter(x_DD_plot, y_DD_plot, 5, 'r');
scatter(x_Wa_new(idx_Wa), y_Wa_new(idx_Wa), 50, 'g', 'filled');
title('After Coordinate Replacement');
legend('Felix\_Wa', 'Felix\_DD', 'Modified Common IDs');
axis equal; grid on;
xlim([-3 3]); ylim([-4.5 4.5]);

% After filtering
subplot(1, 3, 3);
[x_Wa_filt, y_Wa_filt] = latlon2xy([Felix_Wa_filtered.lat], [Felix_Wa_filtered.lon]);

scatter(x_DD_plot, y_DD_plot, 5, 'r'); hold on;
scatter(x_Wa_filt, y_Wa_filt, 50, 'g', 'filled');
title('After Filtering (Only Common IDs)');
legend('Felix\_DD', 'Filtered Felix\_Wa');
axis equal; grid on;
xlim([-3 3]); ylim([-4.5 4.5]);

% Save the filtered dataset if needed
save('/Users/mczhang/Documents/GitHub/FM4/02-data/A_wavelarge5_filtered.mat', 'Felix_Wa_filtered');

% Display final information
fprintf('\nFinal filtered dataset contains only events with these IDs:\n');
disp(commonIDs');
% clc; clear; close all;
% 
% % Load the data
% load('/Users/mczhang/Documents/GitHub/FM4/02-data/Kaiwen_eventData_website.mat')
% Felix_DD = eventData; clear Felix;
% 
% load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_wavelarge5.mat');
% Felix_Wa = Felix; clear Felix;
% 
% % Find common IDs
% commonIDs = intersect([Felix_DD.eventID], [Felix_Wa.ID]);
% 
% % If no common IDs found
% if isempty(commonIDs)
%     error('No common IDs found between Felix_DD and Felix_Wa');
% end
% 
% % Get indices of common events in both datasets
% [~, idx_DD] = ismember(commonIDs, [Felix_DD.eventID]);
% [~, idx_Wa] = ismember(commonIDs, [Felix_Wa.ID]);
% 
% % Convert lat/lon to x/y coordinates for both datasets
% [x_DD, y_DD] = latlon2xy([Felix_DD.lat], [Felix_DD.lon]);
% [x_Wa, y_Wa] = latlon2xy([Felix_Wa.lat], [Felix_Wa.lon]);
% 
% % Extract coordinates for common events
% common_x_DD = x_DD(idx_DD);
% common_y_DD = y_DD(idx_DD);
% common_z_DD = [Felix_DD(idx_DD).depth];
% 
% common_x_Wa = x_Wa(idx_Wa);
% common_y_Wa = y_Wa(idx_Wa);
% common_z_Wa = [Felix_Wa(idx_Wa).depth];
% 
% % Calculate differences in x, y, z coordinates
% dx = common_x_DD - common_x_Wa;
% dy = common_y_DD - common_y_Wa;
% dz = common_z_DD - common_z_Wa;
% 
% % Calculate Euclidean distances
% distances_xy = sqrt(dx.^2 + dy.^2);  % Horizontal distance
% distances_xyz = sqrt(dx.^2 + dy.^2 + dz.^2);  % 3D distance
% 
% % Create figure with histograms
% figure('Position', [100, 100, 1200, 800]);
% 
% % Histogram of x differences
% subplot(2, 3, 1);
% histogram(dx, 200, 'FaceColor', 'b');
% xlabel('Δx (km)');
% ylabel('Count');
% title('X Coordinate Differences');
% grid on;
% 
% % Histogram of y differences
% subplot(2, 3, 2);
% histogram(dy, 200, 'FaceColor', 'g');
% xlabel('Δy (km)');
% ylabel('Count');
% title('Y Coordinate Differences');
% grid on;
% 
% % Histogram of z differences
% subplot(2, 3, 3);
% histogram(dz, 200, 'FaceColor', 'r');
% xlabel('Δz (km)');
% ylabel('Count');
% title('Depth Differences');
% grid on;
% 
% % Histogram of horizontal distances
% subplot(2, 3, 4);
% histogram(distances_xy, 200, 'FaceColor', 'c');
% xlabel('Horizontal Distance (km)');
% ylabel('Count');
% title('Horizontal Distances (XY plane)');
% grid on;
% 
% % Histogram of 3D distances
% subplot(2, 3, 5);
% histogram(distances_xyz, 200, 'FaceColor', 'm');
% xlabel('3D Distance (km)');
% ylabel('Count');
% title('3D Distances (XYZ)');
% grid on;
% 
% % Scatter plot of x vs y differences
% subplot(2, 3, 6);
% scatter(dx, dy, 1, 'filled');
% xlabel('Δx (km)');
% ylabel('Δy (km)');
% title('X vs Y Differences');
% grid on;
% axis equal;
% 
% % Add some statistics to the figure
% annotation('textbox', [0.7, 0.15, 0.2, 0.1], 'String', ...
%     sprintf(['Statistics for %d common events:\n' ...
%              'Mean Δx: %.2f ± %.2f km\n' ...
%              'Mean Δy: %.2f ± %.2f km\n' ...
%              'Mean Δz: %.2f ± %.2f km\n' ...
%              'Mean XY dist: %.2f ± %.2f km\n' ...
%              'Mean 3D dist: %.2f ± %.2f km'], ...
%              length(commonIDs), ...
%              mean(dx), std(dx), ...
%              mean(dy), std(dy), ...
%              mean(dz), std(dz), ...
%              mean(distances_xy), std(distances_xy), ...
%              mean(distances_xyz), std(distances_xyz)), ...
%     'FitBoxToText', 'on', 'BackgroundColor', 'white');
% 
% % Plot the locations of common events
% figure;
% subplot(1, 2, 1);
% scatter(x_DD, y_DD, 5, 'k'); hold on;
% scatter(common_x_DD, common_y_DD,1, 'r', 'filled');
% title('Felix\_DD Events (Common in Red)');
% axis equal;
% xlim([-3 3]);
% ylim([-4.5 4.5]);
% grid on;
% 
% subplot(1, 2, 2);
% scatter(x_Wa, y_Wa, 5, 'k'); hold on;
% scatter(common_x_Wa, common_y_Wa, 1, 'r', 'filled');
% title('Felix\_Wa Events (Common in Red)');
% axis equal;
% xlim([-3 3]);
% ylim([-4.5 4.5]);
% grid on;