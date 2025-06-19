clc;clear;close all;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/G_FM/G_FM_SKHASH_realEle.mat')

 event=event1;clear event1;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/G_FM/G_FM_SKHASH_realEle_manyoutputs.mat')
% event1=event;
for i=1:length(event)
    ind=find([event1.id]==event(i).id);
    if ~isempty(ind) && ind(1)<=5600
    %[rotangle,theta,phi]=kagan([event(i).avmech(1),event(i).avmech(2),event(i).avmech(3)],[event1(ind(1)).avmech(1),event1(ind(1)).avmech(2),event1(ind(1)).avmech(3)]);
    %[rotangle,theta,phi]=kagan([event(i).avmech(1),event(i).avmech(2),event(i).avmech(3)],[event1(ind(1)).strike,event1(ind(1)).dip,event1(ind(1)).rake]);
    [rotangle,theta,phi]=kagan([event(i).strike,event(i).dip,event(i).rake],[event1(ind(1)).strike,event1(ind(1)).dip,event1(ind(1)).rake]);
    
    event(i).kg=rotangle;
    event(i).qual=event1(ind(1)).mechqual;
    clear rotangle;
    end   
end
%event([event.mechqual]=='C' | [event.mechqual]=='D')=[];
%event([event.qual]=='C' | [event.qual]=='D' | [event.qual]=='B')=[];
%event([event.qual]=='C' | [event.qual]=='D')=[];
%event([event.max_azimgap]>100)=[];
figure;set(gcf,'Position',[1000         720         770         518]);
hist(real([event.kg]),120);
title(['Kagan angle difference with 1D and 3D, mean: ' num2str(mean(real([event.kg])),2)]);
%xlim([0 120]);
xlabel('angle')
set(gca,'FontSize',20);
data = real([event.kg]); % Assuming event.kg is your data
% Your existing data and histogram code
% Your existing data and histogram code


% Define the value you are interested in
value = 20;

% Calculate the percentage of data less than or equal to the value
percentage = sum(data <= value) / length(data) * 100;

% Display the percentage
disp(['The percentage of data less than or equal to ', num2str(value), ' is ', num2str(percentage), '%.']);

% Calculate the 50th, 60th, 70th, 80th, and 90th percentiles
p50 = prctile(data, 50);
p60 = prctile(data, 60);
p70 = prctile(data, 70);
p80 = prctile(data, 80);
p90 = prctile(data, 90);

% Display the 50th, 60th, 70th, 80th, and 90th percentile values on the plot
hold on; % Keep the histogram on the plot
yLimits = ylim; % Get current y-axis limits
% 
% % Plot a vertical line at the 50th percentile
% plot([p50 p50], yLimits, 'c--', 'LineWidth', 2); 
% text(p50, yLimits(2) * 0.7, ['50th Percentile: ', num2str(p50)], 'FontSize', 14, 'Color', 'cyan', 'HorizontalAlignment', 'right');
% 
% % Plot a vertical line at the 60th percentile
% plot([p60 p60], yLimits, 'm--', 'LineWidth', 2); 
% text(p60, yLimits(2) * 0.75, ['60th Percentile: ', num2str(p60)], 'FontSize', 14, 'Color', 'magenta', 'HorizontalAlignment', 'right');
% 
% % Plot a vertical line at the 70th percentile
% plot([p70 p70], yLimits, 'b--', 'LineWidth', 2); 
% text(p70, yLimits(2) * 0.8, ['70th Percentile: ', num2str(p70)], 'FontSize', 14, 'Color', 'blue', 'HorizontalAlignment', 'right');
% 
% % Plot a vertical line at the 80th percentile
% plot([p80 p80], yLimits, 'g--', 'LineWidth', 2); 
% text(p80, yLimits(2) * 0.85, ['80th Percentile: ', num2str(p80)], 'FontSize', 14, 'Color', 'green', 'HorizontalAlignment', 'right');
% 
% % Plot a vertical line at the 90th percentile
% plot([p90 p90], yLimits, 'r--', 'LineWidth', 2); 
% text(p90, yLimits(2) * 0.9, ['90th Percentile: ', num2str(p90)], 'FontSize', 14, 'Color', 'red', 'HorizontalAlignment', 'right');

% Plot a vertical line at the specified value
plot([value value], yLimits, 'k--', 'LineWidth', 2); 
text(value, yLimits(2) * 0.95, [num2str(value), ' (', num2str(percentage), '%)'], 'FontSize', 14, 'Color', 'black', 'HorizontalAlignment', 'right');

hold off;

