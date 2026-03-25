clear;close all;
fields={'AS1','AS2','CC1','EC1','EC2','EC3','ID1'};
groups = {'W1','W2','W3','W4','W5','W6','W7','W8','E1','E2','E3','E4','E5','E6','E7','E8'};
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/A_Kaiwen_phase.mat');
Felix([Felix.Pnum]<5)=[];
Felix([Felix.Snum]<5)=[];
dt=1/200; %Hz
%% Step A: obtain wave of event
% window to do filtering
P.a.sttime=-3;P.a.edtime=7;
% window to cut the CC
P.a.window=[-0.25 1];%based on P picks
P.filt=1; % [3 20]
for kz = 1:length(fields)
    for i = 1:length(Felix)
        if isempty(eval(strcat('Felix(i).DDt_', fields{kz})))
            continue;
        end
        trace = obtain_P_waveforms_new(Felix(i), kz, P.a.sttime, P.a.edtime);
        if isempty(trace)
            continue;
        end
        eval(strcat('Felix(i).W_', fields{kz}, '=trace.dataFilt((P.a.window(1)-P.a.sttime)/dt:(P.a.window(2)-P.a.sttime)/dt-1,P.filt);'));
        clear trace;
    end
end
toc
save(['/Users/mczhang/Documents/GitHub/FM3/02-data/A_',groups{gp},'.mat'], 'Felix');
