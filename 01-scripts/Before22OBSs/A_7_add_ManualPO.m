clear; close all;
path = '/Users/mczhang/Documents/GitHub/FM4/02-data/';
path_old = '/Users/mczhang/Documents/GitHub/FM3/02-data/';
groups = {'W1','W2','S1','E1','E2','E3','E4'};
fields = {'AS1','AS2','CC1','EC1','EC2','EC3','ID1'};

% Initialize an empty struct array to collect all wave data
allWaveData = struct([]);

% First loop: Collect all wave data from different files
for gp = 1:length(groups)
    for kz = 1:length(fields)
        try
            % Load wave data
            load([path_old,'D_man/D_',groups{gp},'_',fields{kz},'_uncon.mat']);
            wave=wave_test;
            % Append to our collection
            if isempty(allWaveData)
                allWaveData = wave;
            else
                allWaveData = [allWaveData, wave];
            end
        catch ME
            warning('Failed to load file for group %s field %s: %s', ...
                   groups{gp}, fields{kz}, ME.message);
        end
    end
end

load('/Users/mczhang/Documents/GitHub/FM3/02-data/A_ID/A_group_ALL.mat');
[C,IA,IB]=intersect([Felix.ID],[allWaveData.ID]);
clear allWaveData;

allWaveData=Felix(IA);
% Now load the Felix_Wa_filtered data once
load([path,'A_All/A_wavelarge5_filtered.mat'], 'Felix_Wa_filtered');

% Get the field names of Felix_Wa_filtered
targetFields = fieldnames(Felix_Wa_filtered);

% Initialize the merged data with the original Felix data
mergedData = Felix_Wa_filtered;

% For each entry in the collected wave data
for i = 1:length(allWaveData)
    % Create a new struct with all fields set to 0
    newEntry = struct();
    for j = 1:length(targetFields)
        newEntry.(targetFields{j}) = 0;
    end
    
    % Copy over matching fields from wave
    waveFields = fieldnames(allWaveData);
    for j = 1:length(waveFields)
        fieldName = waveFields{j};
        if isfield(newEntry, fieldName)
            newEntry.(fieldName) = allWaveData(i).(fieldName);
        end
    end
    
    % Add the new entry to mergedData
    mergedData(end+1) = newEntry;
end

% Update Felix_Wa_filtered with the merged data
Felix_Wa_filtered = mergedData;

% Now Felix_Wa_filtered contains:
% 1. All original data from A_wavelarge5_filtered.mat
% 2. All wave data from the various D_..._uncon.mat files
% With missing fields set to 0