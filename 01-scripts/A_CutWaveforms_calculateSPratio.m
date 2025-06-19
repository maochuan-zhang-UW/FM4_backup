clear; close all;

% Define fields and groups
fields = {'AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1'};
groups = {'W1', 'W2', 'W3', 'W4', 'W5', 'W6', 'W7', 'W8', 'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8'};
gp = 2; % Example: set group index (e.g., 'W2'), adjust as needed

% Load data
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_all/A_Kaiwen_phase2F.mat');
dt = 1/200; % Sampling rate (Hz)

%% Load parameters
parameterFile = 'parameter_script_realtimeVer1_MC_focal';
run(parameterFile); % Loads p into workspace

%% Step A: Define parameters for waveform processing
P.a.sttime = -3; % Start time relative to pick
P.a.edtime = 7;  % End time relative to pick
P.a.window = [-0.25 1]; % Window for P-wave extraction
P.filt = 1; % Filter flag (assumed [4 50] Hz, adjust if needed)

% Precompute time vector for amplitude calculations
t_trace = linspace(P.a.sttime, P.a.edtime, 2001)';

% Preallocate NSP fields
for i = 1:length(Felix)
    for kz = 1:length(fields)
        if ~isfield(Felix(i), ['NSP_',fields{kz}])
            Felix(i).(['NSP_', fields{kz}]) = [0, 0, 0];
        end
    end
end

%% Ensure required functions are accessible to parfor workers
% Add the directory containing A_CutWaveforms_calculateSPratio.m to the path
addpath('/Users/mczhang/Documents/GitHub/FM4');
% Verify the file exists
if ~exist('A_CutWaveforms_calculateSPratio.m', 'file')
    error('Cannot find A_CutWaveforms_calculateSPratio.m. Please ensure it is in the MATLAB path.');
end

%% Initialize parallel pool
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool; % Start parallel pool if not already running
end
% Attach the required file to the pool
addAttachedFiles(pool, '/Users/mczhang/Documents/GitHub/FM4/A_CutWaveforms_calculateSPratio.m');

%% Step B: Process amplitudes and waveforms
tic;
for kz = 1:length(fields)
    parfor i = 1:length(Felix)
        % Initialize trace_Z to ensure it’s defined
        trace_Z = [];
        
        % Skip if DDt is empty
        if isempty(Felix(i).(['DDt_', fields{kz}]))
            continue;
        end

        % Step B1: Amplitude calculations
        if isempty(Felix(i).(['NSP_', fields{kz}])) || all(Felix(i).(['NSP_', fields{kz}]) == [0, 0, 0])
            S_pick = Felix(i).(['DDSt_', fields{kz}]) - Felix(i).(['DDt_', fields{kz}]);
            if S_pick >= 0
                % Get waveforms for Z, N, E components
                try
                    [trace_Z, trace_N, trace_E] = obtain_waveforms_all(Felix(i), kz, P.a.sttime, P.a.edtime, p);
                    if ~isempty(trace_Z) && ~isempty(trace_N) && ~isempty(trace_E)
                        % Define windows
                        noise_win = [-0.7, -0.1];
                        P_win = [-0.05, 0.25];
                        S_win = [-0.1, 0.6];
                        
                        if abs(S_pick) < 0.35
                            shift = (0.35 - abs(S_pick)) / 2;
                            P_win = P_win - shift;
                            S_win = S_win + shift;
                        end
                        
                        % Logical indexing
                        noise_idx = t_trace >= noise_win(1) & t_trace <= noise_win(2);
                        P_idx = t_trace >= P_win(1) & t_trace <= P_win(2);
                        S_idx = t_trace >= (S_pick + S_win(1)) & t_trace <= (S_pick + S_win(2));
                        
                        % Extract waveforms
                        Z_noise = trace_Z.dataFilt(noise_idx, 1);
                        N_noise = trace_N.dataFilt(noise_idx, 1);
                        E_noise = trace_E.dataFilt(noise_idx, 1);
                        Z_P = trace_Z.dataFilt(P_idx, 1);
                        N_S = trace_N.dataFilt(S_idx, 1);
                        E_S = trace_E.dataFilt(S_idx, 1);
                        
                        % Calculate amplitudes
                        amp_noise_Z = max(Z_noise) - min(Z_noise);
                        amp_noise_N = max(N_noise) - min(N_noise);
                        amp_noise_E = max(E_noise) - min(E_noise);
                        amp_P_Z = max(Z_P) - min(Z_P);
                        amp_S_N = max(N_S) - min(N_S);
                        amp_S_E = max(E_S) - min(E_S);
                        
                        noise_amp = sqrt((amp_noise_Z.^2 + amp_noise_N.^2 + amp_noise_E.^2) / 3);
                        P_amp = amp_P_Z;
                        S_amp = sqrt((amp_S_N.^2 + amp_S_E.^2) / 2);
                        
                        % Store amplitudes
                        Felix(i).(['NSP_', fields{kz}]) = [noise_amp, S_amp, P_amp];
                    else
                        Felix(i).(['NSP_', fields{kz}]) = [1, 1, 1];
                    end
                catch ME
                    warning('Error in obtain_waveforms_all for Felix(%d), field %s: %s', i, fields{kz}, ME.message);
                    Felix(i).(['NSP_', fields{kz}]) = [1, 1, 1];
                    trace_Z = []; % Ensure trace_Z remains empty
                end
            else
                Felix(i).(['NSP_', fields{kz}]) = [1, 1, 1];
            end
        end

        % Step B2: P-wave waveform extraction (using trace_Z from Step B1)
        if ~isempty(trace_Z)
            try
                % Extract filtered waveform within P.a.window
                idx_start = round((P.a.window(1) - P.a.sttime) / dt);
                idx_end = round((P.a.window(2) - P.a.sttime) / dt);
                Felix(i).(['W_', fields{kz}]) = trace_Z.dataFilt(idx_start:idx_end-1, P.filt);
            catch ME
                warning('Error extracting waveform for Felix(%d), field %s: %s', i, fields{kz}, ME.message);
            end
        end
    end
    fprintf('Elapsed time after processing field %s: %.2f seconds\n', fields{kz}, toc);
end

%% Save output
save_path = '/Users/mczhang/Documents/GitHub/FM4/02-data/A_all/A_Wave2F_7station.mat';
save(save_path, 'Felix');
fprintf('Data saved to %s\n', save_path);