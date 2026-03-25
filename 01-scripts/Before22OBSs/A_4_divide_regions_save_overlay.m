clc; clear; close all;

% Load the original data to get the centroids
load('/Users/mczhang/Documents/GitHub/FM3/02-data/A_All/Felix_kmean_morethan5.mat');

% Remove elements where SP_All is less than 5
Felix([Felix.SP_All] < 5) = [];
[x, y] = latlon2xy([Felix.lat], [Felix.lon]);
XY = [x', y'];

% Define clustering parameters
numClusters = 7;
initialCentroids = [-1, -1; -1.2, -1.8; 0.4,-3.5;1.4,-2.6;1.2,-1.5;1,0.5;1.1,1.3];
namesCluster = {'W1','W2','S1','E4','E3','E2','E1'};

% Perform k-means clustering
[idx, C] = kmeans(XY, numClusters, 'Start', initialCentroids);
clear Felix; 

% Load new dataset
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/A_wavelarge5_filtered_Man_pick.mat');
newData = Felix_Wa_filtered;

% Convert new data's lat/lon to Cartesian coordinates
[x_new, y_new] = latlon2xy([newData.lat], [newData.lon]);
XY_new = [x_new', y_new'];

% Assign each new point to the nearest centroid
distances = pdist2(XY_new, C);
[~, idx_new] = min(distances, [], 2);

% === Save original clusters ===
for i = 1:numClusters
    clusterData = newData(idx_new == i);
    filename = sprintf('/Users/mczhang/Documents/GitHub/FM4/02-data/A_ID/A_%s.mat', namesCluster{i});
    Felix = clusterData;
    save(filename, 'Felix');
    clear Felix;
    fprintf('Saved cluster %d (%s) to %s\n', i, namesCluster{i}, filename);
end

% === Define boundary regions (0.5 km around Voronoi edges) ===
% === Define boundary regions (0.5 km around Voronoi edges) ===
[V, C_voro] = voronoin(C);

boundaryPairs = {
    {'E1', 'E2'}, {'E2', 'E3'}, {'E3', 'E4'}
};
boundaryName ={'E12','E23','E34'};
for pairIdx = 1:length(boundaryPairs)
    pair = boundaryPairs{pairIdx};
    cluster1_idx = find(strcmp(namesCluster, pair{1}));
    cluster2_idx = find(strcmp(namesCluster, pair{2}));
    
    % Debug output
    fprintf('\nProcessing boundary %s-%s (cluster indices %d-%d)\n', ...
            pair{1}, pair{2}, cluster1_idx, cluster2_idx);
    
    % Get events from these two clusters only
    isInEitherCluster = (idx_new == cluster1_idx) | (idx_new == cluster2_idx);
    XY_pair = XY_new(isInEitherCluster,:);
    pairData = newData(isInEitherCluster);
    
    % Find shared Voronoi vertices
    cell1 = C_voro{cluster1_idx};
    cell2 = C_voro{cluster2_idx};
    commonVertices = intersect(cell1, cell2);
    
    if isempty(commonVertices)
        warning('No Voronoi boundary found between %s and %s', pair{1}, pair{2});
        continue;
    end
    
    ridgePoints = V(commonVertices, :);
    
    % Calculate distances to boundary
    if size(ridgePoints, 1) > 1
        boundaryDistances = inf(size(XY_pair, 1), 1);
        for i = 1:size(ridgePoints, 1)-1
            segmentDistances = distanceToLineSegment(XY_pair, ridgePoints(i,:), ridgePoints(i+1,:));
            boundaryDistances = min(boundaryDistances, segmentDistances);
        end
    else
        boundaryDistances = pdist2(XY_pair, ridgePoints);
    end
    
    isNearBoundary = boundaryDistances < 0.5;
    
    % Debug output
    fprintf('Found %d events near %s-%s boundary\n', sum(isNearBoundary), pair{1}, pair{2});
    
    if sum(isNearBoundary) == 0
        warning('No events found within 0.5 km of %s-%s boundary', pair{1}, pair{2});
    end
    
    % Save boundary data
    boundaryData = pairData(isNearBoundary);
    %boundaryName = sprintf('%s_%s', pair{1}, pair{2});
    filename = sprintf('/Users/mczhang/Documents/GitHub/FM4/02-data/A_ID/A_%s.mat', boundaryName{pairIdx});
    
    if ~isempty(boundaryData)
        Felix = boundaryData;
        save(filename, 'Felix');
        clear Felix;
        fprintf('Saved %d events to %s\n', length(boundaryData), filename);
    else
        fprintf('No events to save for %s\n', boundaryName{pairIdx});
    end
end

% Visualization (unchanged from original)
figure; set(gcf, 'position', [0, 500, 650, 779]);
scatter(XY_new(:,1), XY_new(:,2), 1, idx_new, 'filled');
axis equal;
xlim([-3 3]);
ylim([-4.5 4.5]);

hold on;
plot(C(:,1), C(:,2), 'kx', 'MarkerSize', 15, 'LineWidth', 3);
title('New Data with Original Cluster Boundaries');

[vx, vy] = voronoi(C(:,1), C(:,2));
plot(vx, vy, 'k-');

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

% Helper function to calculate distance to line segment
function dist = distanceToLineSegment(points, v1, v2)
    % Vector from v1 to v2
    v = v2 - v1;
    
    % Vector from v1 to each point
    w = bsxfun(@minus, points, v1);
    
    % Projection parameter
    c1 = sum(w .* repmat(v, size(w,1), 1), 2);
    c2 = sum(v .* v);
    
    % Clamp projection to segment
    b = c1 ./ c2;
    b = max(0, min(1, b));
    
    % Projected points
    projection = bsxfun(@plus, v1, bsxfun(@times, b, v));
    
    % Distances
    dist = sqrt(sum((points - projection).^2, 2));
end
