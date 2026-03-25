clc; clear; close all;

% Load the original data to get the centroids (or save/load the centroids directly)
load('/Users/mczhang/Documents/GitHub/FM3/02-data/A_All/Felix_kmean_morethan5.mat');

% Remove elements where SP_All is less than 5 (if needed for centroid calculation)
Felix([Felix.SP_All] < 5) = [];
[x, y] = latlon2xy([Felix.lat], [Felix.lon]);
XY = [x', y'];

% Define the same initial centroids and parameters as in the original clustering
numClusters = 7;
initialCentroids = [-1, -1; -1.2, -1.8; 0.4,-3.5;1.4,-2.6;1.2,-1.5;1,0.5;1.1,1.3];
namesCluster = {'W1','W2','S1','E4','E3','E2','E1'};

% Perform k-means to get the centroids (or just use initialCentroids if no refinement was done)
[idx, C] = kmeans(XY, numClusters, 'Start', initialCentroids);
clear Felix; 

% Now load your NEW dataset
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/A_wavelarge5_filtered_Man_pick.mat') % Example: Replace with your new data
newData = Felix_Wa_filtered; % Assume new_dataset has the same structure as Felix

% Convert new data's lat/lon to Cartesian coordinates
[x_new, y_new] = latlon2xy([newData.lat], [newData.lon]);
XY_new = [x_new', y_new'];

% Assign each new point to the nearest centroid
distances = pdist2(XY_new, C); % Distance from each new point to each centroid
[~, idx_new] = min(distances, [], 2); % Assign to closest centroid

% === Save each cluster as a separate .mat file ===
for i = 1:numClusters
    % Extract events belonging to the current cluster
    clusterData = newData(idx_new == i);
    
    % Create the filename using namesCluster
    filename = sprintf('/Users/mczhang/Documents/GitHub/FM4/02-data/A_ID/A_%sV2.mat', namesCluster{i});

    Felix=clusterData;
    
    % Save the cluster data
    save(filename, 'Felix');
    clear Felix;
    
    fprintf('Saved cluster %d (%s) to %s\n', i, namesCluster{i}, filename);
end

% (Optional) Visualization code remains the same as before)
figure; set(gcf, 'position', [0, 500, 650, 779]);
scatter(XY_new(:,1), XY_new(:,2), 1, idx_new, 'filled');
axis equal;
xlim([-3 3]);
ylim([-4.5 4.5]);

% Overlay original centroids and Voronoi boundaries
hold on;
plot(C(:,1), C(:,2), 'kx', 'MarkerSize', 15, 'LineWidth', 3);
title('New Data with Original Cluster Boundaries');

% Plot Voronoi boundaries
[vx, vy] = voronoi(C(:,1), C(:,2));
plot(vx, vy, 'k-');

% Add other elements (caldera rim, stations, etc.)
axial_calderaRim;
[calderaRim(:,2), calderaRim(:,1)] = latlon2xy(calderaRim(:,2), calderaRim(:,1));
plot(calderaRim(:,2), calderaRim(:,1), 'k', 'LineWidth', 3);

sta = axial_stationsNewOrder;
sta = sta(1:7); 
for i = 1:length(sta)
    [sta(i).x, sta(i).y] = latlon2xy_no_rotate([sta(i).lat], [sta(i).lon]);
    plot(sta(i).x, sta(i).y, 's', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
    text(sta(i).x, sta(i).y, sta(i).name(3:end));
end

grid on;
set(gca, 'FontSize', 12);
xlabel('X km');
ylabel('Y km');