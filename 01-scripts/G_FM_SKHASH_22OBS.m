% Define parameters
clc;clear;close all;
addpath '/Users/mczhang/Documents/GitHub/FM/01-scripts/matlab_ww'
path='/Users/mczhang/Documents/GitHub/FM4/02-data/';
path_FM= '/Users/mczhang/Documents/GitHub/SKHASH/SKHASH';
fields={'AS1','AS2','CC1','EC1','EC2','EC3','ID1','01A', '02A', '03A', '04A', '05A', '06A', '07A','08A','09A','10A', '11A', '12A', '13A', '14A'};
%% input and output
load('/Users/mczhang/Documents/GitHub/FM4/02-data/F_Cl/Filter_Felix_combined.mat');
outputname=[path,'G_FM/G_2F.mat'];
filenameinFM=[path_FM,'/examples/hash3/IN/Axial_22OBSs_Aquality.dat'];
%%

fid=fopen(filenameinFM,'w+');
Po=Po_Clu;
tic;

eventID=[136         246         660        2192        2472        2753        2773        3082        3090        357  4586        4927        5142];
for i=1:length(unique([eventID]))
    ind_Cl=find([Po.cluster]==eventID(i));
% for i=1:length(unique([Po.cluster]))
%     ind_Cl=find([Po.cluster]==i);
    if length(ind_Cl)>=100;continue;end
    %if sum([Po(ind_Cl).ALL])<7;continue;end
    fprintf(fid,'#\n');
    Time_Cl=mean([Po(ind_Cl).on]);
    Lon_Cl=mean([Po(ind_Cl).lon]);
    Lat_Cl=mean([Po(ind_Cl).lat]);
    Depth_Cl=mean([Po(ind_Cl).depth]);
    fprintf(fid,'%3d   %14.7f   %9.4f   %7.4f   %4.2f\n',i,Time_Cl,Lon_Cl,Lat_Cl,Depth_Cl);
    for j=1:length(ind_Cl)
        fprintf(fid,' %6d     %14.7f  %9.4f    %7.4f    %5.3f\n',Po(ind_Cl(j)).ID,Po(ind_Cl(j)).on,Po(ind_Cl(j)).lon,Po(ind_Cl(j)).lat,Po(ind_Cl(j)).depth);
        for k=1:length(fields) % Noise, Pamp, Samp;
            if eval(strcat('Po(ind_Cl(j)).Po_',fields{k},'>0')) %&& eval(strcat('~isempty(Po(ind_Cl(j)).R',fields{k},')')) %if we don't picks, how can we know the polarity
                str_P=['U'];
            elseif eval(strcat('Po(ind_Cl(j)).Po_',fields{k},'<0'))% && eval(strcat('~isempty(Po(ind_Cl(j)).R',fields{k},')'))
                str_P=['D'];
            else
                continue;
            end
            if length(Po(ind_Cl(j)).(['NSP_',fields{k}]))>3
            eval(strcat('Nos=Po(ind_Cl(j)).NSP_',fields{k},'(1);'));%s noise
            eval(strcat('Nop=Po(ind_Cl(j)).NSP_',fields{k},'(2);'));%p noise
            eval(strcat('Pam=Po(ind_Cl(j)).NSP_',fields{k},'(4);'));%p amp
            eval(strcat('Sam=Po(ind_Cl(j)).NSP_',fields{k},'(3);'));%s amp
            elseif length(Po(ind_Cl(j)).(['NSP_',fields{k}]))>2
            eval(strcat('Nos=Po(ind_Cl(j)).NSP_',fields{k},'(1);'));%s noise
            eval(strcat('Nop=Po(ind_Cl(j)).NSP_',fields{k},'(1);'));%p noise
            eval(strcat('Pam=Po(ind_Cl(j)).NSP_',fields{k},'(3);'));%p amp
            eval(strcat('Sam=Po(ind_Cl(j)).NSP_',fields{k},'(2);'));%s amp
            else
            eval(strcat('Nos=1;'));%s noise
            eval(strcat('Nop=1;'));%p noise
            eval(strcat('Pam=1;'));%p amp
            eval(strcat('Sam=1;'));%s amp
            end
            %
            % [ph2dt1(j).([baseName '_yang_nosS']), ...
            %                           ph2dt1(j).([baseName '_yang_nosP']), ...
            %                           ph2dt1(j).([baseName '_yang_Ssnr']), ...
            %                           ph2dt1(j).([baseName '_yang_Psnr'])];
                                      
            letter = char(64 + j);
            str=strcat(letter,fields{k});
            fprintf(fid,[str,'  %s      %10.4f      %10.4f    %12.4f     %12.4f\n'],str_P,Nop,Nos,Pam,Sam);
            %1 sname(1),scomp(1),snet(1),qns1(p),qns2(s),qpamp,qsamp:
        end
    end
end
fprintf(fid,'*\n');
fclose(fid);
J_Write_run_SKHASH_22OBS(filenameinFM,path_FM);
