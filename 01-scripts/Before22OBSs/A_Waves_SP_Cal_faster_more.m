clear; close all;
fields = {'AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1'};
%load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/A_wavelarge5_NSP.mat');
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_ID/A_W2V2.mat')
dt = 1/200; % Hz

%% Load parameters
parameterFile = 'parameter_script_realtimeVer1_MC_focal';
run(parameterFile); % Loads p into workspace

%% Step A: obtain wave of event
P.a.sttime = -3; P.a.edtime = 7;
P.a.window = [-0.25 1];
P.filt = 1; % [4 50] Hz

% Precompute time vector
t_trace = linspace(P.a.sttime, P.a.edtime, 2001)';

% Preallocate NSP fields
for i = 1:length(Felix)
    for kz = 1:length(fields)
        if ~isfield(Felix(i), ['NSP_', fields{kz}])
            Felix(i).(['NSP_', fields{kz}]) = [0, 0, 0];
        end
    end
end

tic;
for kz = 3:length(fields)
    parfor i = 1:length(Felix)
        if isempty(Felix(i).(['NSP_', fields{kz}])) && ~isempty(Felix(i).(['DDt_', fields{kz}]))
            S_pick = Felix(i).(['DDSt_', fields{kz}]) - Felix(i).(['DDt_', fields{kz}]);
            if S_pick >= 0
                % Single call to get all waveforms
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
                    
                   % noise_amp = sqrt((amp_noise_Z^2 + amp_noise_N^2 + amp_noise_E^2) / 3);
                    P_amp = amp_P_Z;
                    %S_amp = sqrt((amp_S_N^2 + amp_S_E^2) / 2);
                    noise_amp = sqrt((amp_noise_Z.^2 + amp_noise_N.^2 + amp_noise_E.^2) / 3);
                    S_amp = sqrt((amp_S_N.^2 + amp_S_E.^2) / 2);
                    
                    % Store amplitudes
                    Felix(i).(['NSP_', fields{kz}]) = [noise_amp, S_amp, P_amp];
                else
                    Felix(i).(['NSP_', fields{kz}]) = [1, 1, 1];
                end
            else
                Felix(i).(['NSP_', fields{kz}]) = [1, 1, 1];
            end
        end
    end
    fprintf('Elapsed time after processing field %s: %.2f seconds\n', fields{kz}, toc);
end

save('/Users/mczhang/Documents/GitHub/FM4/02-data/A_ID/A_W2V2.mat', 'Felix');