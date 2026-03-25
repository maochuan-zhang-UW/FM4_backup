% Define parameters
clc;clear;close all;
addpath '/Users/mczhang/Documents/GitHub/FM/01-scripts/matlab_ww'
path='/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs';
path_FM= '/Users/mczhang/Documents/GitHub/FM/01-scripts/HASH_Manual_5test';
fields={'AS1','AS2','CC1','EC1','EC2','EC3','ID1'};


%% input and output
%load('/Users/mczhang/Documents/GitHub/FM/02-data/E_Cl_toptenW.mat');
load('/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs/F_Cl/F_Cl_All_MLreplace_samecluster_conf.mat');
outputname=[path,'/G_FM/G_HASH_All_ML_sameClusterasbefore.mat'];
filenameinFM=[path_FM,'/Axial_cluster_2026_ML.dat'];

%%

fid=fopen(filenameinFM,'w+');
Po=Po_Clu;
tic;
for i=1:length(unique([Po.Cluster]))
    ind_Cl=find([Po.Cluster]==i);
    if length(ind_Cl)>=100;continue;end
    %if sum([Po(ind_Cl).ALL])<7;continue;end
    fprintf(fid,'#\n');
    Time_Cl=mean([Po(ind_Cl).on]);
    Lon_Cl=mean([Po(ind_Cl).lon]);
    Lat_Cl=mean([Po(ind_Cl).lat]);
    Depth_Cl=mean([Po(ind_Cl).depth]);
    fprintf(fid,'%3d   %14.7f   %9.4f   %7.4f   %4.2f\n',i+99,Time_Cl,Lon_Cl,Lat_Cl,Depth_Cl);
    for j=1:length(ind_Cl)
        fprintf(fid,' %6d     %14.7f  %9.4f    %7.4f    %5.3f\n',Po(ind_Cl(j)).ID,Po(ind_Cl(j)).on,Po(ind_Cl(j)).lon,Po(ind_Cl(j)).lat,Po(ind_Cl(j)).depth);
        for k=1:length(fields) % Noise, Pamp, Samp;
            if eval(strcat('Po(ind_Cl(j)).Po_',fields{k},'(2)>0')) %&& eval(strcat('~isempty(Po(ind_Cl(j)).R',fields{k},')')) %if we don't picks, how can we know the polarity
                str_P=['U'];
            elseif eval(strcat('Po(ind_Cl(j)).Po_',fields{k},'(2)<0'))% && eval(strcat('~isempty(Po(ind_Cl(j)).R',fields{k},')'))
                str_P=['D'];
            else
                continue;
            end
            if length(Po(ind_Cl(j)).(['NSP_',fields{k}]))>3
            eval(strcat('Nos=Po(ind_Cl(j)).NSP_',fields{k},'(1);'));%s noise
            eval(strcat('Nop=Po(ind_Cl(j)).NSP_',fields{k},'(2);'));%p noise
            eval(strcat('Pam=Po(ind_Cl(j)).NSP_',fields{k},'(4);'));%p amp
            eval(strcat('Sam=Po(ind_Cl(j)).NSP_',fields{k},'(3);'));%s amp
            else
            eval(strcat('Nos=Po(ind_Cl(j)).NSP_',fields{k},'(1);'));%s noise
            eval(strcat('Nop=Po(ind_Cl(j)).NSP_',fields{k},'(1);'));%p noise
            eval(strcat('Pam=Po(ind_Cl(j)).NSP_',fields{k},'(3);'));%p amp
            eval(strcat('Sam=Po(ind_Cl(j)).NSP_',fields{k},'(2);'));%s amp
            end
            %
            % [ph2dt1(j).([baseName '_yang_nosS']), ...
            %                           ph2dt1(j).([baseName '_yang_nosP']), ...
            %                           ph2dt1(j).([baseName '_yang_Ssnr']), ...
            %                           ph2dt1(j).([baseName '_yang_Psnr'])];
                                      
            letter = char(64 + j);
            str=strcat(fields{k},letter);
            fprintf(fid,[str,'  %s      %10.4f      %10.4f    %12.4f     %12.4f\n'],str_P,Nop,Nos,Pam,Sam);
            %1 sname(1),scomp(1),snet(1),qns1(p),qns2(s),qpamp,qsamp:
        end
    end
end
fprintf(fid,'*\n');
fclose(fid);
J4_Write_felix_run_HASH3_New(filenameinFM,path_FM);

path_FM= '/Users/mczhang/Documents/GitHub/FM/01-scripts/HASH_Manual_5test';
%cd /Users/mczhang/Documents/GitHub/FM;
filename1=[path_FM,'/hashout1.dat'];
filename2=[path_FM,'/hashout2.dat'];
filename3=[path_FM,'/hashout3.dat'];
[event1] = read_hd3_output1(filename1);
nevent = length(event1);
[event2] = read_hd3_output2(filename2);
nevent2 = length(event2);
[event3] = read_hd3_output3(filename3);
nevent3 = length(event3);
for i=1:length(event1)
    % Change coordinate systems
    [event1(i).avfnorm,event1(i).avslip] =...
        fp2fnorm(event1(i).avmech(:,1),event1(i).avmech(:,2),event1(i).avmech(:,3));
    event1(i).b_axis=cross(event1(i).avfnorm,event1(i).avslip);
    event1(i).t_axis=(event1(i).avfnorm+event1(i).avslip)/sqrt(2);
    event1(i).p_axis=(event1(i).avfnorm-event1(i).avslip)/sqrt(2);
    event1(i).p_l=asin(norm(event1(i).p_axis(3))/norm(event1(i).p_axis))/pi*180;
    event1(i).t_l=asin(norm(event1(i).t_axis(3))/norm(event1(i).t_axis))/pi*180;
    event1(i).b_l=asin(norm(event1(i).b_axis(3))/norm(event1(i).b_axis))/pi*180;
    event1(i).check=sin(event1(i).b_l).^2+sin(event1(i).t_l).^2+sin(event1(i).p_l).^2;
    if event1(i).p_l>= 52 && event1(i).t_l<= 35
        event1(i).faultType = 'N';event1(i).color='b';event1(i).color2=[0,0,1];
    elseif event1(i).p_l<=35 && event1(i).t_l>= 52
        event1(i).faultType = 'R';event1(i).color='r';event1(i).color2=[1,0,0];
    elseif event1(i).p_l<= 40 && event1(i).b_l>= 45 && event1(i).t_l<=40
        event1(i).faultType = 'S';event1(i).color='g';event1(i).color2=[0,1,0];
    elseif event1(i).p_l<= 20 && event1(i).t_l>=40 && event1(i).t_l<= 52
        % olique reverse
        event1(i).faultType = 'R';event1(i).color='r';event1(i).color2=[1,0,0];
    elseif event1(i).p_l>= 40 && event1(i).p_l<= 52 && event1(i).t_l<=20
        % Oblique normal
        event1(i).faultType = 'N';event1(i).color='b';event1(i).color2=[0,0,1];
    else
        event1(i).faultType = 'U';event1(i).color='k';event1(i).color2=[0,0,0];
    end
end
clear event;
%load('/Users/mczhang/Documents/GitHub/FM4/02-data/F_Cl/F_Cl_All.mat');
load('/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs/F_Cl/F_Cl_All_MLreplace_samecluster_conf.mat');
event=event1;
for i=1:length(event1)
    ind=find([Po_Clu.Cluster]==event(i).id-99);
    event(i).lat=mean([Po_Clu(ind).lat]);
    event(i).lon=mean([Po_Clu(ind).lon]);
    event(i).depth=mean([Po_Clu(ind).depth]);
%    event(i).Mw=nanmean([Po_Clu(ind).Mw]);
end
event1=event;
save(outputname, 'event1', 'event2','event3');
toc;
load handel;
sound(y, Fs);