%% Plot and pick the waveform
clc; clear; close all;
dt = 1/200;
P_x = [-0.25 0.25];
P_x1 = P_x(1):dt:P_x(2)-dt;
path = '/Users/mczhang/Documents/GitHub/FM4/02-data/';
fields = {'AS1','AS2','CC1','EC1','EC2','EC3','ID1'};
groups = {'All'};

for gp = 1%:length(groups)
    for kz = 1:length(fields)
        load([path,'C_SVD/C_',fields{kz},'.mat']);
        load([path,'A_ID/A_Wave2F_7station.mat']);
        [commonRows, index5, index2] = intersect([Felix.ID], SVD_result(:,1));
        SVD_result = SVD_result(index2,:);
        [S,I] = sort(abs(SVD_result(:,2)),"descend");
        SVD_result = SVD_result(I,:);
        
        SVD_man = SVD_result(1:min(250, size(SVD_result,1)),:);
        [commonRows, index52, index2] = intersect([Felix.ID], SVD_man(:,1));
        wave = Felix(index52);
        
        non_zero_count = 0;
        for j = 1:length(wave)
            figure(1); clf; set(gcf,'position',[600,500,800,400]);
            eval(strcat('a = wave(j).W_',fields{kz},';'));
            a = a(1:100);
            h(1) = plot(P_x1, normalize(a,'range',[-1 1]), 'b', 'LineWidth', 3);
            hold on
            title(['No:' num2str(wave(j).ID) ' ' fields{kz} ]);
            hold on
            h(3) = plot([0,0], [-0.5,0.5], 'r', 'LineWidth', 2);
            h(2) = plot([-0.0375,-0.0375], [-1,1], 'g', 'LineWidth', 2);
            h(4) = plot([0.0375,0.0375], [-1,1], 'g', 'LineWidth', 2);
            h(5) = plot([-0.1083,-0.1083], [-1,1], 'y', 'LineWidth', 4);
            h(6) = plot([-0.1791,-0.1791], [-1,1], 'k', 'LineWidth', 5);
            h(5) = plot([0.1083,0.1083], [-1,1], 'y', 'LineWidth', 4);
            h(6) = plot([0.1791,0.1791], [-1,1], 'k', 'LineWidth', 5);
            
            [x,y] = ginput(1);
            if x > 0.0375 && x < 0.1083; P = 1;
            elseif x >= 0.1083 && x < 0.1791; P = 2;
            elseif x >= 0.1791 && x < 0.25; P = 3;
            elseif x <= -0.0375 && x > -0.1083; P = -1;
            elseif x <= -0.1083 && x >= -0.1791; P = -2;
            elseif x <= -0.1791 && x > -0.25; P = -3;
            elseif x > -0.0375 && x < 0.0375; P = 0;
            elseif x <= -0.25; P = -4;
            elseif x >= 0.25; P = 4; 
            end
            
            eval(strcat('wave(j).',fields{kz},'_Po = P;'));
            clear a;
            
            % Count non-zero P values
            if P ~= 0
                non_zero_count = non_zero_count + 1;
                if non_zero_count >= 5
                    break; % Exit the loop once we have 5 non-zero P values
                end
            end
        end
        
        % Save the results for this field
        save([path,'D_man/D_',fields{kz},'.mat'], 'wave');
        
        % If we have enough non-zero P values, move to next field
        if non_zero_count >= 5
            continue;
        end
    end
end