% Load combined data
clc;clear;
fields={'AS1','AS2','CC1','EC1','EC2','EC3','ID1','01A', '02A', '03A', '04A', '05A', '06A', '07A','08A','09A','10A', '11A', '12A', '13A', '14A'};

load('/Users/mczhang/Documents/GitHub/FM4/02-data/E_Po/Felix_combined.mat');
load('/Users/mczhang/Documents/GitHub/FM4/02-data/A_all/Kaiwen_eventData_2F.mat');
%[I,IA,IB]=intersect([eventData.eventID],[Felix.ID]);
for i=1:length(Felix_combined)
    Felix_combined(i).depth=eventData(i).depth;
end


% Step 1: Filter based on depth and PoALL
depth_cond = arrayfun(@(x) x.depth > 0.5 && x.depth < 2.5, Felix_combined);
poall_cond = arrayfun(@(x) x.PoALL > 8, Felix_combined);
selected_idx = find(depth_cond & poall_cond);

Filter_Felix_combined = Felix_combined(selected_idx);


for j = 1:length(fields)
    eval(strcat("Felix_combined = rmfield(Felix_combined, 'W_" ,fields{j}, "');"));
end
% Step 4: Add cluster field (1 to N)
for i = 1:length(Filter_Felix_combined)
    Filter_Felix_combined(i).cluster = i;
end

% Save the filtered result
Po_Clu=Filter_Felix_combined;
save('/Users/mczhang/Documents/GitHub/FM4/02-data/F_Cl/Filter_Felix_combined.mat', 'Po_Clu');
