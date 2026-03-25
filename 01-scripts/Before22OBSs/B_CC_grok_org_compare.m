% Compare two filteredResultsMatrix matrices
clc; clear; close all;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/B_CC/B_W1_AS1_grok.mat');
filteredResultsMatrix_grk=filteredResultsMatrix;clear filteredResultsMatrix;
load('/Users/mczhang/Documents/GitHub/FM4/02-data/B_CC/B_W1_AS1.mat');
filteredResultsMatrix_org=filteredResultsMatrix;clear filteredResultsMatrix;
% Assume matrices are already loaded in workspace

% Step 1: Find common rows based on first two columns
[~, idx_grk, idx_org] = intersect(filteredResultsMatrix_grk(:, 1:2), ...
                                 filteredResultsMatrix_org(:, 1:2), 'rows');

% Step 2: Extract common rows
common_grk = filteredResultsMatrix_grk(idx_grk, :);
common_org = filteredResultsMatrix_org(idx_org, :);

% Step 3: Compute differences in the third column
diff_third_column = common_grk(:, 3) - common_org(:, 3);

% Step 4: Summarize differences
num_common_rows = length(diff_third_column);
mean_diff = mean(diff_third_column);
std_diff = std(diff_third_column);
max_abs_diff = max(abs(diff_third_column));
nonzero_diff_count = sum(abs(diff_third_column) > 1e-6); % Count non-zero differences (accounting for floating-point precision)

% Display summary
fprintf('Comparison Summary:\n');
fprintf('Number of common rows: %d\n', num_common_rows);
fprintf('Mean difference (grk - org) in third column: %.6f\n', mean_diff);
fprintf('Standard deviation of differences: %.6f\n', std_diff);
fprintf('Maximum absolute difference: %.6f\n', max_abs_diff);
fprintf('Number of rows with non-zero differences: %d\n', nonzero_diff_count);


% Step 5: Visualize differences
figure;
histogram(diff_third_column, 50, 'EdgeColor', 'k', 'FaceColor', 'b');
title('Histogram of Differences in Third Column (grk - org)');
xlabel('Difference (grk - org)');
ylabel('Frequency');
grid on;

% Step 6: Save results (optional)
save('matrix_comparison_results.mat', 'common_grk', 'common_org', 'diff_third_column');

% Step 7: Display sample of rows with largest differences (if any)
if nonzero_diff_count > 0
    [~, sort_idx] = sort(abs(diff_third_column), 'descend');
    top_n = min(10, nonzero_diff_count); % Show top 10 or fewer
    fprintf('\nTop %d rows with largest absolute differences:\n', top_n);
    fprintf('ID1\t\tID2\t\tgrk_corr\torg_corr\tDifference\n');
    for i = 1:top_n
        idx = sort_idx(i);
        fprintf('%.4f\t%.4f\t%.6f\t%.6f\t%.6f\n', ...
                common_grk(idx, 1), common_grk(idx, 2), ...
                common_grk(idx, 3), common_org(idx, 3), diff_third_column(idx));
    end
end