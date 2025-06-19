% MATLAB script to generate a station correction file similar to north3.statcor.txt
% Output format: <station_code> HHZ OO 0.0000
% Station codes: AS1, AS2, CC1, EC1, EC2, EC3, ID1, 01A, 02A, ..., 14A (21 rows total)
clear;

% MATLAB script to generate a station metadata file similar to scsn.stations.txt
% Output format: <station_code> HHZ <site> <lat> <lon> <elev> <start_date> <end_date> <location>
% Station codes: AS1, AS2, CC1, EC1, EC2, EC3, ID1, 01A, 02A, ..., 14A (21 rows total)

% Define station codes

station_codes = {'01A', '02A', '03A', '04A', '05A', '06A', '07A', ...
                 '08A', '09A', '10A', '11A', '12A', '13A', '14A',...
                 'AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1', ...
                 };

% Define latitude and longitude for each station
latitudes = [45.9336, 45.9338, 45.9547, 45.9496, 45.9397, 45.9361, 45.9257, ...
             46.01933, 46.000625, 45.992922, 45.988743, 45.97439, 45.982938, 45.969067, ...
             45.964147, 45.97084, 45.960338, 45.94959, 45.915367, 45.907817, 45.899965];

longitudes = [-129.9992, -130.0141, -130.0089, -129.9797, -129.9738, -129.9785, -129.9780, ...
              -130.005375, -130.022257, -129.992782, -130.048395, -129.977972, -130.014128, -130.029852, ...
              -130.004503, -130.062393, -129.949987, -130.033977, -130.021773, -129.971555, -130.007573];

% Fixed values for other columns (based on scsn.stations.txt)
channel = 'HHZ';
site = 'BIG CHUCKAWALLA MTNS';
elevation = '   0';
start_date = '1997/09/19';
end_date = '3000/01/01';
location_code = 'OO';

% Output file name
output_file = '/Users/mczhang/Documents/GitHub/FM4/01-scripts/HASH/station.dat';

% Open file for writing
fid = fopen(output_file, 'w');
if fid == -1
    error('Cannot open file %s for writing', output_file);
end

% Write each row in the specified format
for i = 1:length(station_codes)
    % Format latitude and longitude to match the precision in scsn.stations.txt
    % scsn.stations.txt uses 5 decimal places for lat/lon
    lat_str = sprintf('%.5f', latitudes(i));
    lon_str = sprintf('%.5f', longitudes(i));
    
    % Write the line: station_code HHZ site lat lon elev start_date end_date location
    fprintf(fid, '%-4s %-3s  %8s %10s  %-4s %s\n', ...
            ['A' station_codes{i}], channel, lat_str, lon_str, ...
            elevation, location_code);
end

% Close the file
fclose(fid);

disp(['File ', output_file, ' has been generated successfully with 21 rows.']);
%% Define station codes (7 prefixes + 14 numeric codes)
station_codes = {'01A', '02A', '03A', '04A', '05A', '06A', '07A', ...
                 '08A', '09A', '10A', '11A', '12A', '13A', '14A',...
                 'AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1', ...
                 };

% Output file name
output_file = '/Users/mczhang/Documents/GitHub/FM4/01-scripts/HASH/Acor.dat';

% Open file for writing
fid = fopen(output_file, 'w');
if fid == -1
    error('Cannot open file %s for writing', output_file);
end

% Write each station code with HHZ OO 0.0000
for i = 1:length(station_codes)
    fprintf(fid, 'A%s HHZ OO  0.0000\n', station_codes{i});
end

% Close the file
fclose(fid);

disp(['File ', output_file, ' has been generated successfully with 21 rows.']);