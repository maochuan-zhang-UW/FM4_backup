% Data extracted from the image, including elevation
stations = {
    struct('name', 'AX01A', 'lat', 46.01933, 'lon', -130.005375, 'elev', -1588.7),
    struct('name', 'AX01B', 'lat', 46.01935, 'lon', -130.00526, 'elev', -1584.0),
    struct('name', 'AX02A', 'lat', 46.000625, 'lon', -130.029257, 'elev', -1554.7),
    struct('name', 'AX02B', 'lat', 46.00111, 'lon', -130.02923, 'elev', -1554.0),
    struct('name', 'AX03A', 'lat', 45.992922, 'lon', -129.992782, 'elev', -1486.7),
    struct('name', 'AX03B', 'lat', 45.9925, 'lon', -129.99397, 'elev', -1479.0),
    struct('name', 'AX04A', 'lat', 45.988743, 'lon', -130.048395, 'elev', -1535.7),
    struct('name', 'AX04B', 'lat', 45.98834, 'lon', -130.04896, 'elev', -1531.0),
    struct('name', 'AX05A', 'lat', 45.97439, 'lon', -129.977972, 'elev', -1515.7),
    struct('name', 'AX05B', 'lat', 45.97327, 'lon', -129.97728, 'elev', -1512.0),
    struct('name', 'AX06A', 'lat', 45.982938, 'lon', -130.014128, 'elev', -1579.7),
    struct('name', 'AX06B', 'lat', 45.98178, 'lon', -130.01423, 'elev', -1574.0),
    struct('name', 'AX07A', 'lat', 45.969067, 'lon', -130.029852, 'elev', -1580.7),
    struct('name', 'AX07B', 'lat', 45.9695, 'lon', -130.02963, 'elev', -1575.0),
    struct('name', 'AX08A', 'lat', 45.964147, 'lon', -130.004503, 'elev', -1542.7),
    struct('name', 'AX08B', 'lat', 45.96336, 'lon', -130.00441, 'elev', -1529.0),
    struct('name', 'AX09A', 'lat', 45.97084, 'lon', -130.062393, 'elev', -1485.7),
    struct('name', 'AX09B', 'lat', 45.97194, 'lon', -130.05972, 'elev', -1475.0),
    struct('name', 'AX10A', 'lat', 45.960338, 'lon', -129.949987, 'elev', -1568.7),
    struct('name', 'AX10B', 'lat', 45.95867, 'lon', -129.94982, 'elev', -1563.0),
    struct('name', 'AX11A', 'lat', 45.94959, 'lon', -130.039377, 'elev', -1430.7),
    struct('name', 'AX11B', 'lat', 45.94851, 'lon', -130.03801, 'elev', -1426.0),
    struct('name', 'AX12A', 'lat', 45.915367, 'lon', -130.021777, 'elev', -1515.7),
    struct('name', 'AX12B', 'lat', 45.91473, 'lon', -130.02173, 'elev', -1508.0),
    struct('name', 'AX13A', 'lat', 45.907817, 'lon', -129.971555, 'elev', -1564.7),
    struct('name', 'AX13B', 'lat', 45.90472, 'lon', -129.9724, 'elev', -1566.0),
    struct('name', 'AX14A', 'lat', 45.899965, 'lon', -130.007573, 'elev', -1587.7),
    struct('name', 'AX14B', 'lat', 45.89913, 'lon', -130.00612, 'elev', -1585.0),
    struct('name', 'AX15A', 'lat', 45.91972, 'lon', -129.93942, 'elev', -1586.0)
};

% Extract latitude, longitude, and elevation using cellfun
latitudes = cellfun(@(x) x.lat, stations);
longitudes = cellfun(@(x) x.lon, stations);
elevations = cellfun(@(x) x.elev, stations);
% Fix for station names: set UniformOutput to false
station_names = cellfun(@(x) x.name, stations, 'UniformOutput', false);

% Create 3D scatter plot
ax=figure;set(gcf,"Position",[   744    50   885   900])
scatter(longitudes, latitudes,  50, 'blue', 'filled');

% Annotate each point with the station name
for i = 1:length(stations)
    text(longitudes(i), latitudes(i), elevations(i), station_names{i}, ...
         'FontSize', 8, 'HorizontalAlignment', 'right');
end

hold on;
axial_calderaRim;
plot(calderaRim(:,1), calderaRim(:,2), 'k', 'LineWidth', 3);
 hold on;
%%
sta = axial_stationsNewOrder;
sta = sta(1:7); 
for i = 1:length(sta)
    plot(sta(i).lon, sta(i).lat, 's', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
    text(sta(i).lon, sta(i).lat, sta(i).name(3:end));
end
% Add labels and title
title('Stations by Latitude, Longitude, and Elevation');
xlabel('Longitude');
ylabel('Latitude');
zlabel('Elevation (m)');
grid on;


% Save the plot
%saveas(gcf, 'station_3d_scatter_plot.matlab.png');

% Read the data from the file
filename = '/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/nll.temp.202209.202308.txt';  % Replace with your actual filename
data = readtable(filename, 'Format', '%f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 1);

% Extract latitude, longitude, and depth
latitude = data.Var8;
longitude = data.Var9;
depth = data.Var10;

hold on;
scatter3(longitude, latitude,-depth, 5, 'b', 'filled');
title('Earthquake Events');
xlabel('Longitude');
ylabel('Latitude');
grid on;
lonLim = [-130.08, -129.92];
latLim = [45.88, 46.03];
ax_handle = gca;
pbaspect(ax_handle, [diff(lonLim) * cosd(mean(latLim)), diff(latLim), 1]);




%saveas(gcf, 'station_3d_scatter_plot.matlab.png');