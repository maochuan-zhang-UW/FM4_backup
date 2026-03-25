% % % % Read the data
% % % clc;clear;%close all;
% % % %data = dlmread('/Users/mczhang/Documents/GitHub/FM4/02-data/Kaiwen_websit.txt');
% % % 
% % % 
% % % filename = '/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/nll.temp.202209.202308.txt';  % Replace with your actual filename
% % % data = readtable(filename, 'Format', '%f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 1);
% % % data = table2array(data);
% % % % Convert year, month, day, hour, minute, second to datenum
% % % datenum_vals = datenum(data(:,2), data(:,3), data(:,4), ...
% % %                       data(:,5), data(:,6), data(:,7));
% % % 
% % % % Create a structure with cell arrays
% % % eventData = struct(...
% % %     'datenum', num2cell(datenum_vals), ...
% % %     'lat', num2cell(data(:,8)), ...
% % %     'lon', num2cell(data(:,9)), ...
% % %     'depth', num2cell(data(:,10)), ...
% % %     'rms', num2cell(data(:,11)), ...
% % %     'nPha', num2cell(data(:,12)), ...
% % %     'eventID', num2cell(data(:,1)));
% % % 
% % % % Save to a .mat file
% % % save('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/Kaiwen_eventData_2F.mat', 'eventData');  % Cell array version
% % % % OR
% % % %save('eventDataArray.mat', 'eventDataArray');  % Numeric array version

% Load the file
load('A_Kaiwen_phase2Flarger5.mat');

% Check if Felix and the 'on' field exist
if isempty(Felix) || ~isfield(Felix, 'on')
    error('Felix is empty or lacks "on" field');
end

% Extract the 'on' timestamps (datenumbers)
on_times = [Felix.on];

% Verify that on_times is numeric (datenumbers)
if ~isnumeric(on_times)
    error('Felix.on does not contain datenumbers (numeric values). Check the data type.');
end

% Convert datenumbers to datetime
%on_times = datetime(on_times, 'ConvertFrom', 'datenum');

% Convert to dates (truncate time to get only the date part)
dates = dateshift(on_times, 'start', 'day');

% Get unique dates
unique_dates = unique(dates);

% Count events per day
event_counts = histcounts(dates, [unique_dates; unique_dates(end) + days(1)]);

% Display results
disp('Date       | Events');
disp('---------------------');
for i = 1:length(unique_dates)
    fprintf('%s | %d\n', datestr(unique_dates(i), 'dd-mmm-yyyy'), event_counts(i));
end

% Plot histogram of daily counts
figure;
bar(unique_dates, event_counts, 'FaceColor', [0.2 0.6 0.8]);
xlabel('Date');
ylabel('Number of Events');
title('Daily Event Count');
grid on;
datetick('x', 'dd-mmm-yyyy', 'keepticks');
xtickformat('dd-MMM');
xtickangle(45);
