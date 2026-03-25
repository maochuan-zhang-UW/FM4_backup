%% Cross-correlation analysis
clc;clear;close all;
path='/Users/mczhang/Documents/GitHub/FM4/02-data/';
fields={'AS1','AS2','CC1','EC1','EC2','EC3','ID1'};
groups={'W1','W2','S1','E1','E2','E3','E4','W12','E12','E23','E34'};%
dt_window = [0.2,0.4];% short window 0.2 s and long window 0.4 s
windowstring = {'s','l'};
% the frequency of data
dt=1/200; %Hz
% window to cut the CC
P.a.window=[-0.25 1];%based on P picks
% save the CC pairs larger than the threshold below;
P.b.ccthreshold = 0.7;
% 0.1s before the P pick when doing CC
P.b.lagahead = 0.1;
% CC threshold used for the matrix
P.c.ccthreshold=0.7;
% timelag of CC that need to remove,%second
P.c.timelag=0.04;
CC_threshold =P.c.ccthreshold;
TL_threshold =P.c.timelag;

for gp=1:length(groups)
    tic;
    load([path,'/A_ID/A_',groups{gp},'.mat']);
    for kz=1:length(fields)
        for wd=1:length(dt_window)
            % Initialize MAT
            MAT = zeros(1, floor((dt_window(wd)+0.1)/dt));
            n = 1;
            % Populate MAT
            for i = 1:length(Felix)
                fieldName = strcat('Felix(i).W_',fields{kz});
                if ~isempty(eval(fieldName))
                    signal = eval(fieldName);
                    %MAT(n,:)=signal((window(1)-0.1-sttime)/dt:(window(1)-sttime)/dt+dt_window(kz)/dt-1);
                    MAT(n,:)=signal(floor((abs(P.a.window(1))-P.b.lagahead)/dt):floor((abs(P.a.window(1)))/dt)+dt_window(wd)/dt-1);% [-0.1 to the pick]
                    ID_T(n) = Felix(i).ID;
                    n = n + 1;
                end
            end
            % Check if parallel pool exists, if not, create one with 12 workers
            if isempty(gcp('nocreate'))
                parpool(8);
            end
            % Define your signal matrix
            signalMatrix = MAT;

            % Get the length of your time series and the number of time series
            timeSeriesLength = length(signalMatrix(1,:));
            numberOfTimeSeries = length(signalMatrix(:,1));

            % Pre-allocate a cell array to hold the results
            % numberOfTimeSeries*(numberOfTimeSeries-1)/2 is the number of pairs in the worst case
            resultsCell = cell(numberOfTimeSeries*(numberOfTimeSeries-1)/2, 1);

            % Define correlation threshold
            correlationThreshold = P.b.ccthreshold;

            % Initialize a counter for the results
            resultCounter = 1;

            % Calculate cross-correlations for each pair
            for i =1:(numberOfTimeSeries-1)
                tempResults = cell(numberOfTimeSeries - i, 1);
                if isempty(gcp('nocreate'))
                    parpool(8);
                end
                parfor j = (i+1):numberOfTimeSeries
                    % Compute cross-correlation
                    [correlation, lags] = crosscorr(signalMatrix(i,:), signalMatrix(j,:));
                    % Find maximum and minimum correlation
                    [maxCorrelation, maxCorrelationIndex] = max(correlation);
                    [minCorrelation, minCorrelationIndex] = min(correlation);
                    % If max or min correlation is larger than the threshold, store the pair
                    if abs(maxCorrelation) > correlationThreshold || abs(minCorrelation) > correlationThreshold
                        tempResults{j-i} = struct('series1', i, 'series2', j,'series3', ID_T(i), 'series4', ID_T(j), 'maxCorr', maxCorrelation, 'maxLag', lags(maxCorrelationIndex)*dt, 'minCorr', minCorrelation, 'minLag', lags(minCorrelationIndex)*dt);
                    end
                end
                % Remove empty cells from tempResults
                tempResults = tempResults(~cellfun(@isempty, tempResults));
                % Add tempResults to resultsCell
                resultsCell(resultCounter:(resultCounter+length(tempResults)-1)) = tempResults;
                resultCounter = resultCounter + length(tempResults);

            end
            toc;
            % Remove empty cells
            resultsCell = resultsCell(~cellfun(@isempty, resultsCell));
            % Convert resultsCell to a matrix
            resultsMatrix = zeros(length(resultsCell), 8);
            for k = 1:length(resultsCell)
                resultsMatrix(k, :) = [resultsCell{k}.series1, resultsCell{k}.series2,resultsCell{k}.series3, resultsCell{k}.series4, resultsCell{k}.maxCorr, resultsCell{k}.maxLag, resultsCell{k}.minCorr, resultsCell{k}.minLag];
            end

            % Extract series1, series2, maxCorr, and minCorr from resultsMatrix
            series1 = resultsMatrix(:, 1);
            series2 = resultsMatrix(:, 2);
            maxCorr = resultsMatrix(:, 5);
            minCorr = resultsMatrix(:, 7);

            % Compute the maximum value between absolute maxCorr and absolute minCorr, keeping the original sign
            resultsMatrix(:,9) = arrayfun(@(x, y) sign(x * (abs(x) >= abs(y)) + y * (abs(y) > abs(x))) * abs(abs(x) - abs(y)), maxCorr, minCorr);
            % Create a logical index of where the values in column 5 are greater than the absolute values in column 7
            idx = resultsMatrix(:,5) > abs(resultsMatrix(:,7));

            % For those rows where the values in column 5 are greater, save the values from column 5 in column 10 and the corresponding values from column 6 in column 11
            resultsMatrix(idx,10) = resultsMatrix(idx,5);
            resultsMatrix(idx,11) = resultsMatrix(idx,6);

            % For those rows where the absolute values in column 7 are greater or equal, save the original values from column 7 in column 10 and the corresponding values from column 8 in column 11
            resultsMatrix(~idx,10) = resultsMatrix(~idx,7);
            resultsMatrix(~idx,11) = resultsMatrix(~idx,8);
            resultsMatrix(abs(resultsMatrix(:,11)) >= TL_threshold,:)=[];
            resultsMatrix(:,1:2)=[];
            resultsMatrix(:,3:6)=[];
            eval(strcat('resultsMatrix','_',windowstring{wd},'=resultsMatrix;'));
        end
        [commonRows, index5, index2] = intersect(resultsMatrix_l(:,1:2), resultsMatrix_s(:,1:2), 'rows');
        resultsMatrix_cb=resultsMatrix_l(index5,:);
        resultsMatrix_cb(:,6)=resultsMatrix_s(index2,4);
        resultsMatrix_cb(:,7)=resultsMatrix_s(index2,5);
        index3= abs(resultsMatrix_cb(:,7)-resultsMatrix_cb(:,5))>0.01;
        resultsMatrix_cb(index3,:)=[];
        %resultsMatrix_cb_bad=resultsMatrix_cb(index3,:);
        
        filteredResultsMatrix = resultsMatrix_cb;
        filteredResultsMatrix(:,4:7)=[];
        save([path,'B_CC/B_',groups{gp},'_',fields{kz},'.mat'], 'filteredResultsMatrix');
       toc
    end
end