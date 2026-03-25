clc;clear;
% Add the path to obspy-matlab
addpath(genpath('/Users/mczhang/Documents/MATLAB/seizmo'));
rehash toolboxcache;
seizmoverbose(true); 
a=readheader('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/output_test.mseed');
% Read the MiniSEED file
stream = readseizmo('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/output.mseed');

% Loop through traces in the stream
for i = 1:length(stream)
    trace = stream(i);
    % Display information
    fprintf('Station: %s | Start Time: %s | Samples: %d\n', ...
        trace.station, datestr(trace.starttime), trace.npts);

    % Plot the waveform
    figure;
    plot(trace.data);
    title(['Seismic Trace from ', trace.station]);
    xlabel('Sample Number');
    ylabel('Amplitude');
end
