clc; clear;

% === Define input station data ===
stations = {
    'AS1',  45.93360, 129.99921;
    'AS2',  45.93380, 130.01410;
    'CC1',  45.95470, 130.00890;
    'EC1',  45.94960, 129.97971;
    'EC2',  45.93970, 129.97380;
    'EC3',  45.93610, 129.97850;
    'ID1',  45.92570, 129.97800;
    'A01',  46.01933, 130.00537;
    'A02',  46.00062, 130.02226;
    'A03',  45.99292, 129.99278;
    'A04',  45.98874, 130.04839;
    'A05',  45.97439, 129.97797;
    'A06',  45.98294, 130.01413;
    'A07',  45.96907, 130.02985;
    'A08',  45.96415, 130.00450;
    'A09',  45.97084, 130.06239;
    'A10',  45.96034, 129.95000;
    'A11',  45.94959, 130.03398;
    'A12',  45.91537, 130.02177;
    'A13',  45.90782, 129.97156;
    'A14',  45.89997, 130.00757;
};

%% === Write to CSV ===
fid = fopen('/Users/mczhang/Documents/GitHub/SKHASH/SKHASH2/examples/smile/IN/stations.csv', 'w');
fprintf(fid, 'network,station,location,channel,sta_correction\n');

for i = 1:size(stations,1)
    name = stations{i,1};
    lat  = stations{i,2};
    lon  = -stations{i,3}; % west is negative
    fprintf(fid, 'OO,%sA,%.5f %.5f,EHZ,0\n', name, lat, lon);
end

fclose(fid);
disp('stations.csv created with one line per station.');
