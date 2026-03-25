clear;clc;
filename = '/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/hypoDD.temp.202209.202308.pha';
% Open the file
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file');
end

% Initialize variables
Felix = struct();
eventCount = 0;
line = fgetl(fid);

while ischar(line)
    % Check if line is an event header (starts with #)
    if startsWith(line, '#')
        eventCount = eventCount + 1;
        
        % Parse event header
        % Format: # YYYY MM DD HH MM SS.SS lat lon depth ? ? ? ? ID
        data = sscanf(line, '# %d %d %d %d %d %f %f %f %f %f %f %f %f %d');
        
        % Convert time to MATLAB datenum
        eventTime = datenum([data(1:5)', floor(data(6))]) + ...
                    (data(6) - floor(data(6)))/(24*60*60);
        
        % Fill basic fields
        Felix(eventCount).ID = data(14);
        Felix(eventCount).on = eventTime;
        Felix(eventCount).lat = data(7);
        Felix(eventCount).lon = data(8);
        Felix(eventCount).depth = data(9);
        
        % Initialize DDt and DDSt fields to zero
        stations = {'AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1','01A', '02A', '03A', '04A', '05A', '06A', '07A','08A', '09A', '10A','11A', '12A', '13A', '14A', '15A'};
        for i = 1:length(stations)
            Felix(eventCount).(['DDt_' stations{i}]) = 0;
            Felix(eventCount).(['DDSt_' stations{i}]) = 0;
        end
        
        % Process phase arrivals for this event
        nextLine = fgetl(fid);
        while ischar(nextLine) && ~startsWith(nextLine, '#')
            if ~isempty(nextLine)
                % Parse phase arrival
                phase = nextLine(end);      
                %if eventCount < 23686
                    arrivalTime = str2num(nextLine(10:22));
                    
                    station = nextLine(5:7);
                %else
                %    arrivalTime = str2num(nextLine(8:22));
                %    station = nextLine(5:7);
                %end
                
                % Assign to appropriate field based on phase
                if phase == 'P'
                    Felix(eventCount).(['DDt_' station]) = arrivalTime;
                elseif phase == 'S'
                    Felix(eventCount).(['DDSt_' station]) = arrivalTime;
                end
            end
            nextLine = fgetl(fid);
        end
        
        % Count non-zero P and S phases
        Pcount = 0;
        Scount = 0;
        for i = 1:length(stations)
            if Felix(eventCount).(['DDt_' stations{i}]) ~= 0
                Pcount = Pcount + 1;
            end
            if Felix(eventCount).(['DDSt_' stations{i}]) ~= 0
                Scount = Scount + 1;
            end
        end
        Felix(eventCount).Pnum = Pcount;
        Felix(eventCount).Snum = Scount;
        
        line = nextLine;
    else
        line = fgetl(fid);
    end
end

% Close the file
fclose(fid);
save('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/A_Kaiwen_phase.mat','Felix');