% Load the .mat file
clc;clear;

%load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/A_Kaiwen_phase.mat')
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/K_W_041525.mat')

% Convert date numbers to formatted date strings
date_strings = datestr(data_K(:,1), 'dd-mm-yyyy HH:MM:SS');

% Create a new table with date strings instead of date numbers
output_table = table(date_strings, data_K(:,2), data_K(:,3), data_K(:,4), ...
    'VariableNames', {'Date', 'Latitude', 'Longitude', 'Depth'});

% Define output text file name
output_file = 'Seismic_Catalog.txt';

% Write table to text file with formatted headers and data
fid = fopen(output_file, 'w');
% Write header
fprintf(fid, '%-20s %-12s %-12s %-8s\n', 'Date', 'Latitude', 'Longitude', 'Depth');
% Write data
for i = 1:height(output_table)
    fprintf(fid, '%-20s %-12.4f %-12.4f %-8.2f\n', ...
        output_table.Date{i}, output_table.Latitude(i), output_table.Longitude(i), output_table.Depth(i));
end
fclose(fid);

disp(['Text file "' output_file '" has been created.']);