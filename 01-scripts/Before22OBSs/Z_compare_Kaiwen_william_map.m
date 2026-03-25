clc; clear; close all

figure; set(gcf,"Position",[1000, 326, 1177,1142])

for p = 1:2
    namestr = {'Kaiwen', 'William'};
    load('/Users/mczhang/Documents/GitHub/FM4/02-data/K_W_041525.mat');

    if p == 1
        hypo_71 = data_K;
    else
        hypo_71 = data_W;
    end
    clear data;
    

    % Date range filtering
    star = datenum(2022,1,1);
    endd = datenum(2024,1,1);
    indrm = find(hypo_71(:,1) < star | hypo_71(:,1) > endd);
    hypo_71(indrm,:) = [];
    %hypo_71(hypo_71(:,4)<0.1,:)=[];
%    Shallow_depthID(:,1)=hypo_71(hypo_71(:,4)<0.1,5);
%    Shallow_depthID(:,2)=hypo_71(hypo_71(:,4)<0.1,4);
% % Open a text file for writing
% fid = fopen('ShallowDepthID.txt', 'w');
% 
% % Write the matrix to the text file
% for i = 1:size(Shallow_depthID, 1)
%     fprintf(fid, '%d %.4f\n', Shallow_depthID(i, 1), Shallow_depthID(i, 2));
% end
% 
% % Close the file
% fclose(fid);

    [x,y] = latlon2xy(hypo_71(:,2), -hypo_71(:,3));

    % Define subplot positions
    if p == 1
        subplot_positions = [1, 3, 5];  % Left column
    else
        subplot_positions = [2, 4, 6];  % Right column
    end

    % Plotting in the defined subplot positions
    for i = 1:3
        subplot(3, 2, subplot_positions(i));

        hold on;
        axial_calderaRim;
        [calderaRim(:,2), calderaRim(:,1)] = latlon2xy(calderaRim(:,2), calderaRim(:,1));
        plot(calderaRim(:,2), calderaRim(:,1), 'k', 'LineWidth', 2);
        axis equal;
        xlim([-3 3]);
        ylim([-4.5 4.5]);
        zlim([-2 0]);
        
        ind = find(x < 0);
        scatter3(x(ind), y(ind), -hypo_71(ind,4), 2, hypo_71(ind,1), 'filled');
        
        sta = axial_stationsNewOrder;
        sta = sta(1:7); 
        for j = 6
            [sta(j).x, sta(j).y] = latlon2xy([sta(j).lat], [sta(j).lon]);
            plot(sta(j).x, sta(j).y, 's', 'MarkerEdgeColor', 'g', ...
                 'MarkerFaceColor', 'g', 'MarkerSize', 10);
            text(sta(j).x, sta(j).y, sta(j).name(3:end))
        end

        xlabel('x');
        if i == 1
            ylabel('y');
        else
            zlabel('depth');
        end
        
        colormap("jet");
        cbar = colorbar;
        cbar.Ticks = linspace(min(hypo_71(ind,1)), max(hypo_71(ind,1)), 5);
        date_labels = datestr(cbar.Ticks, 'mm/dd/yyyy');
        cbar.TickLabels = date_labels;

        % Define different views for each subplot
        if i == 1
            view(0, 90);  % 
             str='Top view';
        elseif i == 2
            view(0, 0); 
            str='South view';
        else
            view(90, 0);  
            str ='East view';
        end
        
        title([namestr{p} ': ' num2str(length(x(ind))) ' EQs from ' datestr(star, 'yyyy-mm-dd') ' to ' datestr(endd, 'yyyy-mm-dd') ','  str]);

        grid on;
        set(gca, 'Color', [0.5, 0.5, 0.5]);
        hold off;
    end

    %clear;
end
