clear;close all;
fields={'AS1','AS2','CC1','EC1','EC2','EC3','ID1'};
%load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/A_wavelarge5_NSP.mat');
%load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_ID/A_S1.mat');
%load('/Users/mczhang/Documents/GitHub/FM4/02-data/D_man/D_S1_AS1.mat');
load('/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs/D_man/D_S1_AS1.mat')
Felix=wave;
dt=1/200; %Hz
parameterFile = 'parameter_script_realtimeVer1_MC_focal';
run(parameterFile); % Loads p into workspace

%% Step A: obtain wave of event
% window to do filtering
P.a.sttime=-3;P.a.edtime=7;
% window to cut the CC
P.a.window=[-0.25 1];%based on P picks
P.filt=1; % [3 20]
for kz = 2%1:length(fields)
    for i = 1:length(Felix)
        % if ~isempty(eval(strcat('Felix(i).NSP_', fields{kz})))
        %     continue;
        % end        
        if isempty(eval(strcat('Felix(i).DDt_', fields{kz})))
            eval(strcat('Felix(i).NSP_',fields{kz},'(1)=1;')); % Noise
            eval(strcat('Felix(i).NSP_',fields{kz},'(2)=1;')); % Pamp
            eval(strcat('Felix(i).NSP_',fields{kz},'(3)=1;')); % Samp
            continue;
        end
        % trace_Z = obtain_P_waveforms_new(Felix(i), kz, P.a.sttime, P.a.edtime);
        % trace_N = obtain_S_waveforms_N(Felix(i), kz, P.a.sttime, P.a.edtime);
        % trace_E = obtain_S_waveforms_E(Felix(i), kz, P.a.sttime, P.a.edtime)
        [trace_Z, trace_N, trace_E] = obtain_waveforms_all(Felix(i), kz, P.a.sttime, P.a.edtime, p);
        
        if isempty(trace_Z) | isempty(trace_N) | isempty(trace_E)
            eval(strcat('Felix(i).NSP_',fields{kz},'(1)=1;')); % Noise
            eval(strcat('Felix(i).NSP_',fields{kz},'(2)=1;')); % Pamp
            eval(strcat('Felix(i).NSP_',fields{kz},'(3)=1;')); % Samp
            continue;
        end
        %% calculate the noise, P amplitude and S amplitude
        % Define time windows relative to P and S picks
        P_pick = 0; % P pick is at t=0 in the trace
        eval(strcat('S_pick = Felix(i).DDSt_',fields{kz},' - Felix(i).DDt_',fields{kz},';')); % S pick relative to P pick (in seconds)        
        % Standard windows
        if S_pick<0
            eval(strcat('Felix(i).NSP_',fields{kz},'(1)=1;')); % Noise
            eval(strcat('Felix(i).NSP_',fields{kz},'(2)=1;')); % Pamp
            eval(strcat('Felix(i).NSP_',fields{kz},'(3)=1;')); % Samp
            continue;
        end
        noise_win = [-0.7, -0.1]; % Noise window relative to P pick
        P_win = [-0.05, 0.25]; % P-wave window relative to P pick
        S_win = [-0.1, 0.6]; % S-wave window relative to S pick
        
        % Adjust windows if P and S picks are too close (< 0.35 s)
        if abs(S_pick) < 0.35
            shift = (0.35 - abs(S_pick)) / 2;
            P_win = P_win - shift; % Move P window earlier
            S_win = S_win + shift; % Move S window later
        end
        
        % Time vector: -3 to 7 seconds, 2001 samples
        t_trace = linspace(P.a.sttime, P.a.edtime, 2001)';
        
        % Convert time windows to sample indices
        noise_idx = find(t_trace >= noise_win(1) & t_trace <= noise_win(2));
        P_idx = find(t_trace >= P_win(1) & t_trace <= P_win(2));
        S_idx = find(t_trace >= (S_pick + S_win(1)) & t_trace <= (S_pick + S_win(2)));
        
        % Extract waveforms for each window from dataFilt(:,1)
        Z_noise = trace_Z.dataFilt(noise_idx, 1);
        N_noise = trace_N.dataFilt(noise_idx, 1);
        E_noise = trace_E.dataFilt(noise_idx, 1);
        
        Z_P = trace_Z.dataFilt(P_idx, 1);
        %N_P = trace_N.dataFilt(P_idx, 1);
        %E_P = trace_E.dataFilt(P_idx, 1);
        
        %Z_S = trace_Z.dataFilt(S_idx, 1);
        N_S = trace_N.dataFilt(S_idx, 1);
        E_S = trace_E.dataFilt(S_idx, 1);
        
        % Calculate max-min amplitudes for each window and channel
        amp_noise_Z = max(Z_noise) - min(Z_noise);
        amp_noise_N = max(N_noise) - min(N_noise);
        amp_noise_E = max(E_noise) - min(E_noise);
        
        amp_P_Z = max(Z_P) - min(Z_P);
        % amp_P_N = max(N_P) - min(N_P); % Not used for P amplitude
        % amp_P_E = max(E_P) - min(E_P); % Not used for P amplitude
        
        amp_S_N = max(N_S) - min(N_S);
        amp_S_E = max(E_S) - min(E_S);
        
        % Calculate final amplitudes
        noise_amp = sqrt((amp_noise_Z^2 + amp_noise_N^2 + amp_noise_E^2) / 3); % RMS of three channels
        P_amp = amp_P_Z; % Vertical channel amplitude
        S_amp = sqrt((amp_S_N^2 + amp_S_E^2) / 2); % RMS of N and E channels
        
        % Store amplitudes in Felix structure
        eval(strcat('Felix(i).NSP_', fields{kz}, '(1)=', num2str(noise_amp), ';'));
        eval(strcat('Felix(i).NSP_', fields{kz}, '(2)=', num2str(S_amp), ';'));
        eval(strcat('Felix(i).NSP_', fields{kz}, '(3)=', num2str(P_amp), ';'));
        
        %%
        % Add this plotting section inside the loop after calculating the windows
        % For demonstration, plot for the first event (i=1) and first channel (kz=1, AS1)
        %if i == 1 && kz == 1 % Plot only for the first event and AS1 channel
            figure('Position', [100, 100, 800, 600]);

            % Define colors for the waveforms
            colors = {'b', 'r', 'g'}; % Z: blue, N: red, E: green

            % Create subplots for Z, N, E channels
            for ch = 1:3
                subplot(3, 1, ch);
                hold on;

                % Select the appropriate trace
                if ch == 1
                    trace = trace_Z.dataFilt(:, 1); % Z channel
                    ch_label = 'Z (Vertical)';
                elseif ch == 2
                    trace = trace_N.dataFilt(:, 1); % N channel
                    ch_label = 'N (North-South)';
                else
                    trace = trace_E.dataFilt(:, 1); % E channel
                    ch_label = 'E (East-West)';
                end

                % Plot the waveform
                plot(t_trace, trace, colors{ch}, 'LineWidth', 1.5, 'DisplayName', ch_label);

                % Highlight windows with shaded regions
                % Noise window
                patch([noise_win(1), noise_win(2), noise_win(2), noise_win(1)], ...
                    [min(trace)*1.2, min(trace)*1.2, max(trace)*1.2, max(trace)*1.2], ...
                    'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'Noise Window');

                % P window
                patch([P_win(1), P_win(2), P_win(2), P_win(1)], ...
                    [min(trace)*1.2, min(trace)*1.2, max(trace)*1.2, max(trace)*1.2], ...
                    'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'P Window');

                % S window
                patch([S_pick + S_win(1), S_pick + S_win(2), S_pick + S_win(2), S_pick + S_win(1)], ...
                    [min(trace)*1.2, min(trace)*1.2, max(trace)*1.2, max(trace)*1.2], ...
                    'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', 'S Window');

                % Add P and S pick lines
                xline(0, 'k--', 'P Pick', 'LineWidth', 1, 'HandleVisibility', 'off');
                xline(S_pick, 'm--', 'S Pick', 'LineWidth', 1, 'HandleVisibility', 'off');

                % Customize plot
                xlabel('Time relative to P pick (s)');
                ylabel('Amplitude');
                title(sprintf('Waveform for Event %d, Channel %s (%s)', i, fields{kz}, ch_label));
                grid on;
                legend('show', 'Location', 'northeast');
                hold off;
            end

            % Adjust layout
            sgtitle(sprintf('Three-Component Waveform for Event %d, Channel %s', i, fields{kz}));
            set(gcf, 'Color', 'w');
            %close;
        %end
        clear trace_Z trace_N trace_E;
    end
    fprintf('Elapsed time after processing field %s: %.2f seconds\n', fields{kz}, toc);
end

%save(['/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/A_wavelarge5_NSP.mat'], 'Felix');
