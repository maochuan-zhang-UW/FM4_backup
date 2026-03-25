%% This script is to find any corresponding waveforms, cut based on the on+ ddt
function [traceEvent]=obtain_S_waveforms_E(ph2dt,ind_station,st1,st2)
%% define the parameter
parameterFile = 'parameter_script_realtimeVer1_MC_focal.m';
eval(parameterFile(1:end-2));
stations={'AXAS1','AXAS2','AXCC1','AXEC1','AXEC2','AXEC3','AXID1'};
channels={'EHE',  'EHE',   'HHE',  'EHE',   'HHE',  'EHE',  'EHE'};
% define the station to plotå
stations=stations(ind_station);
channels=channels(ind_station);
% Get the waveforms
d=ph2dt.on;
fileData = [ p.dir.data '/' datestr(d, 'yyyy') '/' ...
    datestr(d, 'mm') '/' datestr(d, 'yyyy-mm-dd-HH-00-00') '.mat'];
% Define the file path

% Check if the file exists
if exist(fileData, 'file') == 2
    trace1 = load(fileData);
    trace1=trace1.trace;clear trace;% confilct with a function
else
    trace1 =[];
end

% Continue with the rest of your code


p.focal.tlimEvent=[st1 st2];
if ~isempty(trace1)
    % Filter Data and demean [4 50]Hz
    eval(strcat('phaseT=ph2dt.DDt_',stations{1}(3:end),';'));
    %eval(strcat('phaseT=0;'));% MZ change it to get the time based on the on time.
    traceEvent = subset_trace(trace1,ph2dt.on+(phaseT+p.focal.tlimEvent)/86400,stations, channels,[],0);
    clear trace1;
    if length(traceEvent)==0 || length(traceEvent.data)<700
        traceEvent=[];
    else
        for i=1:length(traceEvent)
            for j=1:length(p.filt) % 2 is [4  50] Hz as waldehauser
                try
                    traceEvent(i).dataFilt(:,j) = trace_filter(traceEvent(i).data,p.filt(j),traceEvent(i).sampleRate);
                catch exception
                    keyboard;
                end
                traceEvent(i).dataFilt(:,j) = traceEvent(i).dataFilt(:,j) - mean(traceEvent(i).dataFilt(:,j));
            end
        end
    end
else
    traceEvent=[];
end
end