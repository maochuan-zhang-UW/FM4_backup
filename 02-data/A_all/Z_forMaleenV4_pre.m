% File name
clc;clear;
fname = 'hypo71.dat';

% Open file
fid = fopen(fname,'r');

% Skip header line
fgetl(fid);

% Read the data
C = textscan(fid, ...
    '%8s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f %f', ...
    'Delimiter',' ', 'MultipleDelimsAsOne',true);

fclose(fid);

% Parse date and time
yyyymmdd = C{1};     % cell array of strings
HHMM = C{2};         % hours and minutes
SSS  = C{3};         % seconds (with decimals)

n = numel(yyyymmdd);
S = struct([]);

for i = 1:n
    % Convert date string
    datestr_raw = yyyymmdd{i};
    yyyy = str2double(datestr_raw(1:4));
    mm   = str2double(datestr_raw(5:6));
    dd   = str2double(datestr_raw(7:8));

    % Convert HHMM + seconds
    hh = floor(HHMM(i)/100);
    min = mod(HHMM(i),100);
    sec = SSS(i);

    % MATLAB datenum
    t = datenum(yyyy,mm,dd,hh,min,sec);

    % Lat: deg + min/60
    lat = C{4}(i) + C{5}(i)/60;

    % Lon: deg + min/60
    lon = C{6}(i) + C{7}(i)/60;

    % Depth
    depth = C{8}(i);

    % Mw
    Mw = C{9}(i);

    % ID
    ID = C{16}(i);

    % Save into struct
    S(i).ID    = ID;
    S(i).time  = t;     % datenum now
    S(i).lat   = lat;
    S(i).lon   = lon;
    S(i).depth = depth;
    S(i).Mw    = Mw;
end

% --- Extract arrays ---
lat  = [S.lat]';
lon  = [S.lon]';   % already negative (°W convention)
side = {S.side}';

% --- Masks ---
mask_west = strcmp(side,'w');
mask_east = strcmp(side,'e');

% --- Plot ---
figure; set(gcf,'Position',[400 400 800 600]);

% West (blue)
scatter(lon(mask_west), lat(mask_west), 25, 'b', 'filled'); hold on;

% East (red)
scatter(lon(mask_east), lat(mask_east), 25, 'r', 'filled');

xlabel('Longitude (°W)');
ylabel('Latitude (°N)');
title('Earthquake Locations: West vs East');
legend('West region','East region','Location','best');
grid on;

% Flip longitude axis for °W
set(gca,'XDir','reverse');
