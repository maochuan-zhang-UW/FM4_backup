% Load both files
load('/Users/mczhang/Documents/GitHub/FM4/02-data/E_Po/E_Po_15OBS.mat');  % First structure
Felix1 = Felix;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/E_Po/E_Po.mat');        % Second structure
Felix2 = Felix;

% Make sure they have the same number of events
if length(Felix1) ~= length(Felix2)
    error('The number of events in the two files do not match.');
end

% Initialize combined structure
Felix_combined = Felix1;

% Loop through all events
for i = 1:length(Felix1)
    if Felix1(i).ID ~= Felix2(i).ID
        error('Mismatch in event ID at index %d: %d ≠ %d', i, Felix1(i).ID, Felix2(i).ID);
    end

    fields2 = fieldnames(Felix2(i));
    
    for j = 1:length(fields2)
        field = fields2{j};

        if strcmp(field, 'PoALL')
            % Sum PoALL from both
            Felix_combined(i).PoALL = Felix1(i).PoALL + Felix2(i).PoALL;
        elseif ~isfield(Felix1(i), field)
            % Copy only fields not in Felix1
            Felix_combined(i).(field) = Felix2(i).(field);
        end
    end
end

% Save the result
save('/Users/mczhang/Documents/GitHub/FM4/02-data/E_Po/Felix_combined.mat', 'Felix_combined');
