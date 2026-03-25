% Define parameters
clc;clear;close all;
addpath '/Users/mczhang/Documents/GitHub/FM/01-scripts/matlab_ww'
path='/Users/mczhang/Documents/GitHub/FM3/02-data/';
path_FM= '/Users/mczhang/Documents/GitHub/SKHASH/SKHASH2';
fields={'AS1','AS2','CC1','EC1','EC2','EC3','ID1','A01', 'A02', 'A03', 'A04', 'A05', 'A06', 'A07',...
    'A08', 'A09', 'A10', 'A11', 'A12', 'A13', 'A14'};
%% input and output
load('/Users/mczhang/Documents/GitHub/FM4/02-data/F_Cl/F_Cl_SKHASH3D.mat');
Po=Po_Clu;clear Po_Clu;
tic;
   n=1;m=1;
%ID=[5          13          63         103         113         257         260         292];
ID=unique([Po.cluster]);%[ 257         260         292];
%for i=1:length(unique([Po.Cluster])) 
for i=1:length(unique(ID)) 
    ind_Cl=find([Po.cluster]==ID(i));
    if length(ind_Cl)>=100;continue;end
    %if sum([Po(ind_Cl).ALL])<7;continue;end
    timenum(m)=mean([Po(ind_Cl).on]);
    longitude(m)=mean([Po(ind_Cl).lon]);
    latitude(m)=mean([Po(ind_Cl).lat]);
    depth(m)=mean([Po(ind_Cl).depth]);
    horz_uncertain_km(m)=0.3;
    vert_uncertain_km(m)=0.2;
    mag(m)=0;
    event_id2(m)=ID(i);
    m=m+1;
    for j=1:length(ind_Cl)
        for k=1:length(fields) % Noise, Pamp, Samp;
            if mean(Po(ind_Cl(j)).(['NSP_',fields{k}]))==1
                continue;
            end
            if mean(Po(ind_Cl(j)).(['AZ',fields{k}]))==1
                continue;
            end
            event_id(n)=ID(i);
            letter = char(64 + j);
            station{n}=strcat(fields{k},letter);
            network{n}='XX';
            location{n}='--';
            channel{n}='EHZ';
            if eval(strcat('Po(ind_Cl(j)).Po_',fields{k},'>0')) %&& eval(strcat('~isempty(Po(ind_Cl(j)).R',fields{k},')')) %if we don't picks, how can we know the polarity
                p_polarity(n)=1;
                %p_polarity(n)=Po(ind_Cl(j)).(['Po_', fields{k}]);
            elseif eval(strcat('Po(ind_Cl(j)).Po_',fields{k},'<0'))% && eval(strcat('~isempty(Po(ind_Cl(j)).R',fields{k},')'))
                p_polarity(n)=-1;
                %p_polarity(n)=Po(ind_Cl(j)).(['Po_', fields{k}]);
            else
                continue;
            end
            if length(Po(ind_Cl(j)).(['NSP_',fields{k}]))>3
            eval(strcat('Nos=Po(ind_Cl(j)).NSP_',fields{k},'(1);'));%s noise
            eval(strcat('Nop=Po(ind_Cl(j)).NSP_',fields{k},'(2);'));%p noise
            eval(strcat('Pam=Po(ind_Cl(j)).NSP_',fields{k},'(4);'));%p amp
            eval(strcat('Sam=Po(ind_Cl(j)).NSP_',fields{k},'(3);'));%s amp
            elseif length(Po(ind_Cl(j)).(['NSP_',fields{k}]))==2
                    eval(strcat('Nos=1;'));%s noise
            eval(strcat('Nop=1;'));%p noise
            eval(strcat('Pam=1;'));%p amp
            eval(strcat('Sam=1;'));%s amp

            else
            eval(strcat('Nos=Po(ind_Cl(j)).NSP_',fields{k},'(1);'));%s noise
            eval(strcat('Nop=Po(ind_Cl(j)).NSP_',fields{k},'(1);'));%p noise
            eval(strcat('Pam=Po(ind_Cl(j)).NSP_',fields{k},'(3);'));%p amp
            eval(strcat('Sam=Po(ind_Cl(j)).NSP_',fields{k},'(2);'));%s amp
            end 
            %sp_ratio(n)=log10(Sam/Pam);
            sp_ratio(n)=Sam/Pam;
            a=Po(ind_Cl(j)).(['TOA',fields{k}]);
            if ~isempty(a)
             % eval(strcat('takeoff(n)=Po(ind_Cl(j)).event.TOA',fields{k},';'));
             % eval(strcat('azimuth(n)=Po(ind_Cl(j)).event.AZ',fields{k},';'));
             eval(strcat('takeoff(n)=Po(ind_Cl(j)).TOA',fields{k},';'));
             eval(strcat('azimuth(n)=Po(ind_Cl(j)).AZ',fields{k},';'));
            else 
                takeoff(n)=0;
                azimuth(n)=0;
            end

            % if isfield(Po(ind_Cl(j)), ['TOA', fields{k}]) 
            %   eval(strcat('takeoff(n)=Po(ind_Cl(j)).TOA',fields{k},';'));
            %   eval(strcat('azimuth(n)=Po(ind_Cl(j)).AZ',fields{k},';'));
            % end
            takeoff_uncertainty(n)=10;
            azimuth_uncertainty(n)=5;
            n=n+1;

        end
    end
end

% Remove entries where both takeoff and azimuth are 0
valid_idx = ~(takeoff == 0 & azimuth == 0);

% Filter all arrays accordingly
event_id = event_id(valid_idx);
station = station(valid_idx);
network = network(valid_idx);
location = location(valid_idx);
channel = channel(valid_idx);
sp_ratio = sp_ratio(valid_idx);
p_polarity = p_polarity(valid_idx);
takeoff = takeoff(valid_idx);
takeoff_uncertainty = takeoff_uncertainty(valid_idx);
azimuth = azimuth(valid_idx);
azimuth_uncertainty = azimuth_uncertainty(valid_idx);


% Create a table
T = table(event_id', station', network', location', channel', sp_ratio', takeoff', takeoff_uncertainty', azimuth', azimuth_uncertainty', ...
    'VariableNames', {'event_id', 'station', 'network', 'location', 'channel', 'sp_ratio', 'takeoff', 'takeoff_uncertainty', 'azimuth', 'azimuth_uncertainty'});
% Write the table to a CSV file
writetable(T, '/Users/mczhang/Documents/GitHub/SKHASH/SKHASH2/examples/smile/IN/amp.csv');

% Create a table
T1 = table(event_id', station', network', location', channel', p_polarity', takeoff', takeoff_uncertainty', azimuth', azimuth_uncertainty', ...
    'VariableNames', {'event_id', 'station', 'network', 'location', 'channel', 'p_polarity', 'takeoff', 'takeoff_uncertainty', 'azimuth', 'azimuth_uncertainty'});
% Write the table to a CSV file
writetable(T1, '/Users/mczhang/Documents/GitHub/SKHASH/SKHASH2/examples/smile/IN/pol.csv');

clear event_id;
event_id=event_id2;
    time= datetime(timenum, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');


% Create a table
T2 = table( string(time)', latitude', longitude', depth', horz_uncertain_km', vert_uncertain_km', mag',event_id', ...
    'VariableNames', {'time', 'latitude', 'longitude', 'depth', 'horz_uncertain_km', 'vert_uncertain_km', 'mag','event_id'});

% Write the table to a CSV file
writetable(T2, '/Users/mczhang/Documents/GitHub/SKHASH/SKHASH2/examples/smile/IN/eq_catalog.csv');

cd /Users/mczhang/Documents/GitHub/SKHASH/SKHASH2/
!rm /Users/mczhang/Documents/GitHub/SKHASH/SKHASH2/examples/smile/OUT/out.csv
command = '/opt/miniconda3/envs/SKHASH/bin/python SKHASH.py ./examples/smile/control_file.txt';
[status, result] = system(command);
disp(result);
readSKHASH3D_skhash2('/Users/mczhang/Documents/GitHub/FM4/02-data/G_FM/G_3D.mat')
%HA_SKHASH_read_result;
c