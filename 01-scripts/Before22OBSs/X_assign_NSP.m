%help me to understand the code below:
clear;close all;
path='/Users/mczhang/Documents/GitHub/FM4/02-data/';
path_old='/Users/mczhang/Documents/GitHub/FM3/02-data/';
fields={'AS1','AS2','CC1','EC1','EC2','EC3','ID1'};
groups = {'E1','E2','E3','E4','S1','W1','W2'};
nz=1;
resultman=[];
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/A_wavelarge5_NSP.mat');
Felix_all=Felix; clear Felix;
for gp=1:length(groups)
    load([path,'/A_ID/A_',groups{gp},'.mat']);
    for i = 1:length(Felix)
        for kz = 1:length(fields)
            if ~isfield(Felix(i), ['NSP_', fields{kz}])
                Felix(i).(['NSP_', fields{kz}]) = [0, 0, 0];
            end
        end
    end
    [C,IA,IB]=intersect([Felix.ID],[Felix_all.ID]);
    for i=1:length(C)
        for kz=1:length(fields)
            Felix(IA(i)).(['NSP_', fields{kz}])=Felix_all(IB(i)).(['NSP_', fields{kz}]);
        end
    end
    save(['/Users/mczhang/Documents/GitHub/FM4/02-data/A_ID/A_' groups{gp} '.mat'], 'Felix');
    clear Felix;
end
