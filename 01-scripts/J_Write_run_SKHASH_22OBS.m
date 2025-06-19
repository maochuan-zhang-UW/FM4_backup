function J_Write_run_SKHASH_22OBS(filename,path)
fid = fopen(filename,'r');
Nline=nan;
while 1
    tline = fgetl(fid);
    if strcmp(tline,'*');B(id)=Nline;break;end
    if strcmp(tline,'#')		% then this is an earthquake line
        if ~isnan(Nline)
            B(id)=Nline;
        end
        tline = fgetl(fid);
        Nline=0;
        id=str2num(tline(1:4));
        %if id==199;disp('prepare to pause');pause(5);end
    % elseif strcmp(tline(1:2),'AX') % MZ change for testing new hash
    %     Nline=Nline+1;
    elseif ismember(tline(1), {'A', 'B'})
         Nline=Nline+1;
    end
end
fclose(fid);
%% 
phaFile=[path '/examples/hash3/IN/north2.phase.txt'];
ampFile=[path,'/examples/hash3/IN/north3.amp.txt'];
if exist(phaFile) == 2
    delete(phaFile);
end

if exist(ampFile) == 2
    delete(ampFile);
end

%%
%filename=['Axial_cluster.dat'];
fid = fopen(filename,'r');
fid2=fopen([path '/examples/hash3/IN/north2.phase.txt'],'w+');
fid3=fopen([path,'/examples/hash3/IN/north3.amp.txt'],'w+');
id=nan;
while 1
    tline = fgetl(fid);
    if strcmp(tline,'*')
        fprintf(fid2,['                                                        ' ...
            '%16s\n'],id);
        break;
    end  % this is the end of the file
    if strcmp(tline,'#')		% then this is an earthquake line
        if ~isnan(id)
            fprintf(fid2,['                                                        ' ...
                '%16s\n'],id);
        end
        tline = fgetl(fid);
        id=tline(1:4);
        %if strcmp(tline(1:3),num2stop); break; end
        Amp_num=B(str2num(id));
        fprintf(fid3,'%d         %d\n',str2num(id), Amp_num);
        t= datetime(str2num(tline(7:20)),'ConvertFrom','datenum');
        yr=year(t);mo=month(t);da=day(t);hr=hour(t);mn=minute(t);sec=second(t);
        ilat = floor(abs(str2num(tline(36:42))));
        mlat = 60*(abs(str2num(tline(36:42))) - ilat);
        cns = 'N'; %if sign(location.lat)==-1; cns='S'; end;
        ilon = floor(abs(str2num(tline(24:32))));
        mlon = 60*(abs(str2num(tline(24:32))) - ilon);
        cew = 'W'; %if sign(location.lon)==-1; cew='W'; end;
        dep = str2num(tline(46:49));
        if dep<0.5;dep=0.5;end
        %eh = 0.08;
        %ez = 0.05;
         eh = 0.3;
         ez = 0.2;
        mag = 1;
        fprintf(fid2,...
            '%4i%2i%2i%2i%2i%5.2f%2i %5.2f%3i %5.2f%5.2f     %3i%41s%5.2f  %5.2f   %36s%4.2f%22s\n',...
            yr,mo,da,hr,mn,sec,ilat,mlat,ilon,mlon,dep,B(str2num(id)),'',eh,ez,'',mag,id);
    %     fprintf(fid2,...
    % '%4i%02i%02i%02i %02i %06.2f %10.5f %10.5f %7.2f %49s %4.2f %6.2f %77s %4.2f %22s\n',...
    % yr, mo, da, hr, mn, sec, origin_lat, origin_lon, origin_depth_km, ' ', horz_uncert_km, vert_uncert_km, ' ', magnitude, id);

    elseif ismember(tline(1), {'A'})
        sta=tline(1:3);
        staReal=tline(1:4);
        schan = 'HHZ';
        cpol=tline(7);
        oneset = 'I';
        fprintf(fid2,'%4s %2s  %3s %c %c\n',staReal, ...
            'OO',schan,oneset,cpol);
        Noip=str2num(tline(14:24));
        Nois=str2num(tline(28:40));
        Pamp=str2num(tline(44:56));
        Samp=str2num(tline(60:72));
        fprintf(fid3,'%4s %3s %2s  %5.2f   %5.2f  %10.3f %10.3f %10.3f %10.3f\n',staReal,schan,'OO',Noip,Noip,Noip,Nois,Pamp,Samp);
    end
end
fclose(fid);
fclose(fid2);
fclose(fid3);

%cd /Users/mczhang/Documents/GitHub/FM/01-scripts/HASH_Manual/
%!./hash_driver3 < hash.input 

% data=readmatrix('hashout1.dat');
% FMangle(:,1)=data(:,1);
% FMangle(:,2)=data(:,22);
% FMangle(:,3)=data(:,23);
% FMangle(:,4)=data(:,24);
%addpath /Users/mczhang/Documents/GitHub/Axial-AutoLocate/HASH_code/
%cd ..
%end
end