clc; clear;

% T = readtable('/Users/mczhang/Documents/GitHub/FM4/02-data/Kaiwen_websit.txt');
% T = table2struct(T);
% for i = 1:length(T)
%     data_K(i,1) = datenum([T(i).Var1], [T(i).Var2], [T(i).Var3], [T(i).Var4], [T(i).Var5], [T(i).Var6]);
%     data_K(i,2) = T(i).Var7;
%     data_K(i,3) = -T(i).Var8;
%     data_K(i,4) = T(i).Var9;
%     data_K(i,5) = T(i).Var14;
% end
% data_K2 = data_K; clear data_K T;
% 
% % Read and process the second dataset (hypoDD)
% T = readtable('/Users/mczhang/Documents/GitHub/FM/02-data/Alldata/Axial_Kaiwen_ML_2021.txt');
% T = table2struct(T);
% for i = 1:length(T)
%     data_K(i,1) = datenum([T(i).Var1], [T(i).Var2], [T(i).Var3], [T(i).Var4], [T(i).Var5], [T(i).Var6]);
%     data_K(i,2) = T(i).Var7;
%     data_K(i,3) = -T(i).Var8;
%     data_K(i,4) = T(i).Var9;
%     data_K(i,5) = T(i).Var15;
% end
% data_K = [data_K; data_K2]; clear data_K1 data_K2;
% data_K(data_K(:,1)<datenum(2015,1,22),:)=[];
% 
% % Read and process the third dataset (HYPO71)
% fileHypo71 = '/Users/mczhang/Documents/GitHub/FM4/02-data/hypo71.dat';
% expanded = 1;
% header = 1;
% hypo71 = read_HYPO71(fileHypo71, expanded, header);
% data(:,1) = hypo71.datenum;
% data(:,2) = hypo71.latDeg + hypo71.latMin / 60;
% data(:,3) = hypo71.lonDeg + hypo71.lonMin / 60;
% data(:,4) = hypo71.depth;
% % Your existing code to process the datasets goes here...
% % Omitted for brevity, make sure to include all your data processing steps above

load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_All/K_W_041525.mat');
data=data_W;
% Define the bin edges for the histograms
dday = 30;
bin_edges = min(min(data(:,1)), min(data_K(:,1))):dday:max(max(data(:,1)), max(data_K(:,1)));

% Calculate the histograms
[counts_william, ~] = histcounts(data(:,1), bin_edges);
[counts_kaiwen, ~] = histcounts(data_K(:,1), bin_edges);

% Calculate the ratio of the number of events per day
ratio = counts_kaiwen ./ counts_william;
ratio(isnan(ratio)) = 0;  % Replace NaNs with zeros

% Plotting
% Plotting
figure; set(gcf,'Position',[ 653         595        1158         643]);

% Combined plot for William and Kaiwen
yyaxis left;
b1 = bar(bin_edges(1:end-1), counts_william, 'FaceColor', 'b', 'EdgeColor', 'b', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5); % Adjust color and transparency as needed
hold on;
b2 = bar(bin_edges(1:end-1), counts_kaiwen, 'FaceColor', 'r', 'EdgeColor', 'r', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5); % Adjust color and transparency as needed
ylabel('Event Count');
ylim([0 20000]);

% Ratio plot on the right y-axis
yyaxis right;
p1 = plot(bin_edges(1:end-1), ratio, 'k', 'LineWidth', 2); % Black line for the ratio
ylabel('Ratio (Kaiwen/William)');
ylim([0 10]);
xlim([min(data(:,1)) max(data(:,1))]);
hold on
plot([min(data(:,1)) max(data(:,1))],[1 1],LineWidth=2,Color='g')
hold off;

% Add legend
legend([b1, b2, p1], {'William', 'Kaiwen', 'Ratio (Kaiwen/William)'}, 'Location', 'best');

title(['Eqs comparison, Count every ' num2str(dday) ' days']);
xlabel('Date');

% Adjust x-axis to display dates
specific_dates = {'2016-01-01', '2017-01-01', '2018-01-01', '2019-01-01', '2020-01-01', '2021-01-01', '2022-01-01', '2023-01-01', '2024-01-01','2025-01-01','2026-01-01'};
date_nums = datenum(specific_dates, 'yyyy-mm-dd');
ax = gca;
ax.XTick = date_nums;  % Set the ticks to the specified dates
ax.XTickLabel = datestr(ax.XTick, 'yyyy-mm-dd');  % Set the tick labels to the specified dates format

