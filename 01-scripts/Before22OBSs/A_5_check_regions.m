clear;close all
path='/Users/mczhang/Documents/GitHub/FM4/02-data/';
namesCluster={'W1','W2','S1','E1','E2','E3','E4','W12','E12','E23','E34'};%

namesClusterv2={'R7','R6','R5','R1','R2','R3','R4','R67','R12','R23','R34'};


% Create a figure for all subplots

figure;set(gcf,'Position',[ 1000         499        1038         838])
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/CenterXY.mat');
% Loop over each cluster
for kz=1:length(namesCluster)
    subplot(3,4,kz); % Arrange in a 3x4 grid
    load([path,'A_ID/A_',namesCluster{kz},'.mat']);

    % Filtering
    %Felix([Felix.ID] > 200100019) = [];

    % Convert latitude and longitude to Cartesian coordinates
    [x, y] = latlon2xy([Felix.lat], [Felix.lon]);

    % Combine x and y coordinates into a single matrix
    XY = [x', y'];

    % Scatter plot of the data points colored by cluster index
    % Use a unique color for each subplot
    scatter(XY(:,1), XY(:,2), 1, 'filled', 'CData', rand(1,3)); % Random color

    % Axis settings
    axis equal;
    xlim([-3 3]);
    ylim([-4.5 4.5]);
    hold on;

    % Plot caldera rim
    axial_calderaRim;
    [calderaRim(:,2), calderaRim(:,1)] = latlon2xy(calderaRim(:,2), calderaRim(:,1));
    plot(calderaRim(:,2), calderaRim(:,1), 'k', 'LineWidth', 3);
 hold on;
 [vx, vy] = voronoi(C(:,1), C(:,2));
plot(vx, vy, 'k-');
    % Set title
    set(gca,'FontSize',12)
    title([namesClusterv2{kz}, ', Num: ' num2str(length(XY))]);
   
    % Clear variables to avoid conflicts in the next iteration
    clear Felix B A com;
end

sgtitle('7 regions and 4 overlap regions');
set(gca,'FontSize',12)
hold off;
savemyfigure
