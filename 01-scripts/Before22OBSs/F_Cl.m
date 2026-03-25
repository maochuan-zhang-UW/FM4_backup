clc; clear; close all;

% Constants and paths
BASE_PATH = '/Users/mczhang/Documents/GitHub/';
addpath([BASE_PATH, 'FM/01-scripts/subcode/']);
data_path = [BASE_PATH, 'FM4/02-data/'];
fields = {'AS1','AS2','CC1','EC1','EC2','EC3','ID1'};
groups = {'W','E','All'};
min_Po_num = 4;
min_SP_num = 4;
average_num = 6;
fluc_num = 1;

% Test different 'a' values
%a_values = [0.01, 0.1, 1, 10, 30];
a_values = [  0.2];
%a_values = [ 0.1:0.01: 0.3];
results = cell(length(a_values), 1);

for gp = 3:length(groups) % Process all groups
    % Load and filter data
    %load([BASE_PATH, 'FM/02-data/D_Po_', groups{gp}, '.mat']);
    load('/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs/E_Po/E_All.mat')
    Po = Felix([Felix.PoALL] > min_Po_num & [Felix.SP_All] > min_SP_num);
    numbertoknow = length(Po);

    % Pre-allocate data matrices
    data_location = [vertcat(Po.lat), vertcat(Po.lon), vertcat(Po.depth)];
    [data_location(:,1), data_location(:,2)] = latlon2xy(data_location(:,1), data_location(:,2));

    n_fields = length(fields);
    polarities = zeros(length(Po), n_fields);
    sp_ratios = zeros(length(Po), n_fields);

    % Vectorized data preparation
    for kz = 1:n_fields
        polarities(:,kz) = vertcat(Po.(['Po_' fields{kz}]));
        sp_ratios(:,kz) = vertcat(Po.(['SP_' fields{kz}]));
    end

    % Main data matrix
    data_matrix_base = [data_location, polarities, sp_ratios, vertcat(Po.ID)];

    % Test each 'a' value
    for a_idx = 1:length(a_values)
        tic;
        a = a_values(a_idx);
        b = 1;
        c = 100;

        % Working copy of data
        data_matrix = data_matrix_base;
        data_location_curr = data_location;
        polarities_curr = polarities;
        sp_ratios_curr = sp_ratios;

        % Pre-allocate result storage
        group_matrix = [];
        group_numbers = [];
        loc_values = [];
        Po_values = [];
        Po_nums = [];
        Po_ratio = [];
        Spr_values = [];
        group_sizes = []; % Store group sizes
        group_count = 1;

        % Define the sequence of minimum group sizes
        min_group_sizes = [11, 10, 9, 8, 7, 6, 5, 4, 3];
        current_index = 1;

        while size(data_matrix, 1) > average_num && current_index <= length(min_group_sizes)
            minLeaves = min_group_sizes(current_index);
            maxLeaves = minLeaves + fluc_num; % Allow some flexibility

            % Vectorized distance calculations
            eucD_loc = pdist(data_location_curr, 'seuclidean');
            eucD_spr = pdist(sp_ratios_curr, @custom_distance_SPr);
            eucD_po = pdist(polarities_curr, @custom_distance_Po);
            eucD = a * eucD_loc + b * eucD_spr + c * eucD_po;

            clustTreeEuc = linkage(eucD, 'complete');
            n = size(data_matrix, 1);

            % Efficient cluster size calculation
            numLeaves = ones(2*n-1, 1);
            for i = 1:size(clustTreeEuc, 1)
                numLeaves(n+i) = sum(numLeaves(clustTreeEuc(i,1:2)));
            end

            selectedNodes = find(numLeaves(n+1:end) >= minLeaves & ...
                numLeaves(n+1:end) <= maxLeaves) + n;

            if isempty(selectedNodes)
                current_index = current_index + 1;
                continue;
            end

            % Pre-allocate cluster metrics
            leavesall = cell(length(selectedNodes), 1);
            loc_val = zeros(length(selectedNodes), 1);
            Po_val = zeros(length(selectedNodes), 1);
            Po_nval = zeros(length(selectedNodes), 1);
            Po_raval = zeros(length(selectedNodes), 1);
            Spr_val = zeros(length(selectedNodes), 1);
            DIS = zeros(length(selectedNodes), 1);

            % Process clusters
            for i = 1:length(selectedNodes)
                leaves = get_leaves(clustTreeEuc, selectedNodes(i), n);
                leavesall{i} = leaves;
                group_data = data_matrix(leaves,:);

                % Location metrics
                loc_mean = mean(group_data(:,1:3));
                loc_diff = group_data(:,1:3) - loc_mean;
                loc_val(i) = sqrt(mean(sum(loc_diff.^2, 2)));

                % Polarity metrics
                data_po = group_data(:,4:10);
                [mismatches, non_zeros] = arrayfun(@(col) ...
                    calc_mismatches(data_po(:,col)), 1:size(data_po,2));
                Po_val(i) = sum(mismatches);
                Po_nval(i) = sum(non_zeros);
                Po_raval(i) = Po_val(i)/Po_nval(i);

                % SP ratio metrics
                Spr_val(i) = mean(pdist(group_data(:,11:17), @custom_distance_SPr));

                DIS(i) = a * loc_val(i) + b * Spr_val(i) + c * Po_val(i);
            end

            % Sort clusters by DIS values (ascending)
            [~, sortIdx] = sort(DIS);

            % Initialize containers for selected clusters and their events
            selectedClusters = [];
            assignedEvents = [];

            % Iterate through sorted clusters
            for idx = sortIdx'
                clusterEvents = leavesall{idx};

                % Check for overlap with already assigned events
                if isempty(intersect(clusterEvents, assignedEvents))
                    % Apply criteria
                    if Po_val(idx) == 0 && Spr_val(idx) <= 0.2 && loc_val(idx) <= 0.3
                        selectedClusters = [selectedClusters, idx];
                        assignedEvents = [assignedEvents; clusterEvents(:)];
                    end
                end
            end

            % Process the results
            if isempty(selectedClusters)
                % No clusters selected, move to next smaller minLeaves
                current_index = current_index + 1;
            else
                % Process all selected clusters
                for selIdx = selectedClusters
                    leaveindex = leavesall{selIdx};
                    group_data = data_matrix(leaveindex,:);
                    group_data(:,end+1) = group_count;

                    group_matrix = [group_matrix; group_data];
                    group_numbers = [group_numbers; group_count];
                    loc_values = [loc_values; loc_val(selIdx)];
                    Po_values = [Po_values; Po_val(selIdx)];
                    Po_nums = [Po_nums; Po_nval(selIdx)];
                    Po_ratio = [Po_ratio; Po_raval(selIdx)];
                    Spr_values = [Spr_values; Spr_val(selIdx)];
                    group_sizes = [group_sizes; length(leaveindex)]; % Store group size

                    fprintf('a=%.2f: Grouped %d events into group %d; RMS Loc: %.2f, Po: %d/%d, Spr: %.2f\n', ...
                        a, length(leaveindex), group_count, loc_val(selIdx), Po_val(selIdx), Po_nval(selIdx), Spr_val(selIdx));

                    group_count = group_count + 1;
                end

                % Remove assigned events from the dataset
                data_matrix(assignedEvents, :) = [];
                data_location_curr(assignedEvents, :) = [];
                polarities_curr(assignedEvents, :) = [];
                sp_ratios_curr(assignedEvents, :) = [];
            end
        end

        % Handle remaining events
        if ~isempty(data_matrix)
            data_matrix(:,end+1) = group_count;
            group_matrix = [group_matrix; data_matrix];
            group_sizes = [group_sizes; size(data_matrix, 1)]; % Store remaining group size
        end

        % Calculate total grouped events
        total_grouped_events = size(group_matrix, 1);

        elapsed = toc;
        results{a_idx} = struct('a', a, 'group_matrix', group_matrix, ...
            'group_count', group_count, 'loc_values', loc_values, ...
            'Po_values', Po_values, 'Po_nums', Po_nums, 'Po_ratio', Po_ratio, ...
            'Spr_values', Spr_values, 'group_sizes', group_sizes, ...
            'total_grouped_events', total_grouped_events, 'time', elapsed);

        fprintf('a=%.2f: Grouped %d events into %d clusters in %.2f seconds\n', ...
            a, total_grouped_events, group_count, elapsed);
    end

    % Visualization for all 'a' values
    figure('Position', [200 400 1200 800]);
    sgtitle(['Group: ', groups{gp}]); % Add group name to figure title (MATLAB R2018b+)
    for a_idx = 1:length(a_values)
        r = results{a_idx};
        
        % Subplot 1: Location RMS with group sizes and total events
        subplot(length(a_values), 4, (a_idx-1)*4 + 1);
        bar(r.loc_values);
        title(sprintf('a=%.2f: Loc RMS (Mean: %.2f)', r.a,mean(r.loc_values)));
        % for i = 1:length(r.loc_values)
        %     text(i, r.loc_values(i) + 0.05 * max(r.loc_values), num2str(r.group_sizes(i)), ...
        %         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        % end
        ylim([0, max(r.loc_values) * 1.1]);

        % Subplot 2: Polarity Ratio with total events
        subplot(length(a_values), 4, (a_idx-1)*4 + 2);
        bar(r.Po_ratio);
        title(sprintf('b=100: Po (Group/Tot: %d / %d)',sum(r.group_sizes(1:end-1)), r.total_grouped_events));

        % Subplot 3: S/P Ratio with total events
        subplot(length(a_values), 4, (a_idx-1)*4 + 3);
        bar(r.Spr_values);
        title(sprintf('c=1: S/P Ratio (Mean: %.2f)', mean(r.Spr_values)));

        % Subplot 4: Group Size Distribution
        subplot(length(a_values), 4, (a_idx-1)*4 + 4);
        bar(r.group_sizes(1:end-1),10);
        title(sprintf(' Mean Group Size %.2f', r.a,mean(r.group_sizes(1:end-1))));
        ylabel('Event number');
    end
    % If Po does not have a field named 'Cluster', add it
    if ~isfield(Po, 'Cluster')
        [Po.Cluster] = deal(NaN);  % or some other default value
    end

    % Initialize Po_Clu with the same structure as Po
    Po_Clu = Po;
    Po_Clu = Po_Clu(1:length(group_matrix(:,1)));  % resize to the correct length

    for i=1:length(group_matrix(:,1))
        loc1=find([Po.ID]==group_matrix(i,18));       
        Po_Clu(i)=Po(loc1(1));
        Po_Clu(i).Cluster=group_matrix(i,19);
    end
    %save([path,'F_Cl/F_',groups{gp},'.mat'], 'Po_Clu');
    save('/Users/mczhang/Documents/GitHub/FM4/02-data/F_Cl/F_Cl_All.mat', 'Po_Clu');
    disp([groups{gp},'_group numer:',num2str(length(Po_Clu)),', total:',num2str(numbertoknow),', ratio:',num2str(length(Po_Clu)/numbertoknow)]);
    save(['/Users/mczhang/Documents/GitHub/FM4/02-data/F_Cl/F_Cl_' groups{gp} '_new_Cluster_Result.mat'], "results");
end

% Helper functions
function leaves = get_leaves(clustTree, node, n)
queue = node;
leaves = [];
while ~isempty(queue)
    curr = queue(1);
    queue(1) = [];
    if curr <= n
        leaves = [leaves; curr];
    else
        queue = [queue; clustTree(curr-n, 1:2)'];
    end
end
end

function [mismatches, non_zeros] = calc_mismatches(col)
nz = col(col ~= 0);
if isempty(nz)
    mismatches = 0; non_zeros = 0;
else
    mismatches = sum(nz ~= mode(nz));
    non_zeros = length(nz);
end
end
