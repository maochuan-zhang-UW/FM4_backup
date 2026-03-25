
%% add magnitude to the Po_Clu
clc;clear;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/F_Cl/F_All.mat');
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/Kaiwen_eventData_website.mat')
[C,IA,IB]=intersect([eventData.eventID],[Po_Clu.ID]);
%PMomMagnitude = 2/3*log10(PMoment)-10.7;
for i=1:length(C)
    Po_Clu(IB(i)).Mw=eventData(IA(i)).magnitude;
    %Po_Clu(IB(i)).Mo=10^((3/2) * data(IA(i),2) + 10.7);
end
% 
figure;hist([Po_Clu.Mw],200);
title('The Mw of earthquakes with FMs');
set(gca,'fontSize',16)
save('/Users/mczhang/Documents/GitHub/FM4/02-data/F_Cl/F_All.mat', 'Po_Clu');
% %load('/Users/mczhang/Documents/GitHub/FM3/02-data/G_FM/G_HASH_Final_loc02.mat');

load('/Users/mczhang/Documents/GitHub/FM4/02-data/F_Cl/F_All.mat')
load('/Users/mczhang/Documents/GitHub/FM4/02-data/G_FM/G_HASH_All.mat')
clear event2 event3;
for i=1:length(event1)
    ind=find([Po_Clu.Cluster]==event1(i).id-99);
    for j=1:length(ind)
    Po_Clu(ind(j)).avfnorm=event1(i).avfnorm;
    Po_Clu(ind(j)).avslip= event1(i).avslip;
    Po_Clu(ind(j)).avslip= event1(i).avslip;
    Po_Clu(ind(j)).avmech= event1(i).avmech;
    Po_Clu(ind(j)).mechqual= event1(i).mechqual;
    Po_Clu(ind(j)).color= event1(i).color;
    Po_Clu(ind(j)).color2= event1(i).color2;
    Po_Clu(ind(j)).faultType= event1(i).faultType;
    Po_Clu(ind(j)).max_azimgap= event1(i).max_azimgap;
    Po_Clu(ind(j)).max_takeoff= event1(i).max_takeoff;
    end
end

[B,I]=sort([Po_Clu.on],'ascend');
Po_Clu=Po_Clu(I);



n=1;
for i=1:length(Po_Clu)
    if ~isempty(Po_Clu(i).avfnorm)
        Po_Clu1(n)=Po_Clu(i);
        n=n+1;
    end
end
clear Po_Clu;
Po_Clu=Po_Clu1;
Po_Clu_org=Po_Clu;

Po_Clu([Po_Clu.max_azimgap]>240)=[];
Po_Clu([Po_Clu.max_takeoff]>60)=[];
save('/Users/mczhang/Documents/GitHub/FM4/02-data/G_FM/G_All_finalgap.mat', "Po_Clu",'Po_Clu_org');
% Po_Clu = rmfield(Po_Clu, 'SP_AS1');
% Po_Clu = rmfield(Po_Clu, 'SP_AS2');
% Po_Clu = rmfield(Po_Clu, 'SP_CC1');
% Po_Clu = rmfield(Po_Clu, 'SP_EC1');
% Po_Clu = rmfield(Po_Clu, 'SP_EC2');
% Po_Clu = rmfield(Po_Clu, 'SP_ID1');
% Po_Clu = rmfield(Po_Clu, 'SP_EC3');
% 
% Po_Clu = rmfield(Po_Clu, 'DDt_AS1');
% Po_Clu = rmfield(Po_Clu, 'DDt_AS2');
% Po_Clu = rmfield(Po_Clu, 'DDt_CC1');
% Po_Clu = rmfield(Po_Clu, 'DDt_EC1');
% Po_Clu = rmfield(Po_Clu, 'DDt_EC2');
% Po_Clu = rmfield(Po_Clu, 'DDt_ID1');
% Po_Clu = rmfield(Po_Clu, 'DDt_EC3');
% 
% Po_Clu = rmfield(Po_Clu, 'NSP_AS1');
% Po_Clu = rmfield(Po_Clu, 'NSP_AS2');
% Po_Clu = rmfield(Po_Clu, 'NSP_CC1');
% Po_Clu = rmfield(Po_Clu, 'NSP_EC1');
% Po_Clu = rmfield(Po_Clu, 'NSP_EC2');
% Po_Clu = rmfield(Po_Clu, 'NSP_ID1');
% Po_Clu = rmfield(Po_Clu, 'NSP_EC3');
% 
% Po_Clu = rmfield(Po_Clu, 'Po_AS1');
% Po_Clu = rmfield(Po_Clu, 'Po_AS2');
% Po_Clu = rmfield(Po_Clu, 'Po_CC1');
% Po_Clu = rmfield(Po_Clu, 'Po_EC1');
% Po_Clu = rmfield(Po_Clu, 'Po_EC2');
% Po_Clu = rmfield(Po_Clu, 'Po_ID1');
% Po_Clu = rmfield(Po_Clu, 'Po_EC3');