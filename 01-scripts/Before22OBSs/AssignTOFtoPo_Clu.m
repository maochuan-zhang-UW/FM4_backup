load('AZ_TOF_newCluster.mat');
load('/Users/mczhang/Documents/GitHub/FM3/02-data/F_Cl/F_2015Erp_Final_AZ_TOA_1D23D.mat')
% List of stations
stations = {'AS1','AS2','CC1','EC1','EC2','EC3','ID1'};

% Create a map for fast lookup from event.id to its index
eventIDMap = containers.Map([event.id], 1:numel(event));

for i = 1:numel(Po_Clu)
    if isKey(eventIDMap, Po_Clu(i).ID)
        eidx = eventIDMap(Po_Clu(i).ID);  % find the matching event
        for s = 1:numel(stations)
            st = stations{s};

            % Get existing values (make sure they’re row vectors)
            az_old = Po_Clu(i).(['AZ' st]);
            toa_old = Po_Clu(i).(['TOA' st]);

            % Get the new value from event
            az_new = event(eidx).(['AZ' st]);
            toa_new = event(eidx).(['TOA' st]);

            % Append as the 3rd value
            Po_Clu(i).(['AZ' st]) = [az_old(:); az_new];
            Po_Clu(i).(['TOA' st]) = [toa_old(:); toa_new];
        end
    end
end
