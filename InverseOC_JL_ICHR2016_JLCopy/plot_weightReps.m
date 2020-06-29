% plot_weightoverreps
% plot data
function [h, rsq_mean, rsq_std] = plot_weightReps(h, outputPath, currInstName, feature_full, avgWeightArray_belowThres, ...
    resnromAll_lsqlin_const_minRMSE_array_belowThres, cost_function_names, ...
    segmentInfo, plotGroups, runMode)
ax = [];

% runMode = 'motionOnly_seg'; % motionOnly_seg, startToEnd_seg, motionOnly_rest, startoEnd_rest

% calc mean/std based on the plotgroups and segmentinfo timing
for ind = 1:length(plotGroups)
    % find the corresponding entries
    indsToUse = []; 
    
    if plotGroups{ind} ~= 0
        switch runMode
            case 'startToEnd_seg'
                currStartTime = segmentInfo.timeStart(plotGroups{ind}(1));
                currEndTime = segmentInfo.timeEnd(plotGroups{ind}(end));
                
                [startVal, startInd] = findClosestValue(currStartTime, feature_full.t);
                [endVal, endInd] = findClosestValue(currEndTime, feature_full.t);
                
                indsToUse = startInd:endInd;
                
            case 'motionOnly_seg'
                for ind_seg = 1:length(plotGroups{ind})
                    currStartTime = segmentInfo.timeStart(plotGroups{ind}(ind_seg));
                    currEndTime = segmentInfo.timeEnd(plotGroups{ind}(ind_seg));
                    
                    [startVal, startInd] = findClosestValue(currStartTime, feature_full.t);
                    [endVal, endInd] = findClosestValue(currEndTime, feature_full.t);
                    
                    indsToUse = [indsToUse startInd:endInd];
                end
                
            case 'startToEnd_rest'
                % want to pull out the resting period after the reported
                % fatigue level has been hit
                currStartTime = segmentInfo.timeStart(plotGroups{ind}(end));
                
                switch ind
                    case length(plotGroups)
                        currEndTime = feature_full.t(end);
                        
                    otherwise
                        currEndTime = segmentInfo.timeStart(plotGroups{ind+1}(1));
                end
                
                [startVal, startInd] = findClosestValue(currStartTime, feature_full.t);
                [endVal, endInd] = findClosestValue(currEndTime, feature_full.t);
                
                indsToUse = startInd:endInd;
                
            case 'motionOnly_rest'
                 for ind_seg = 1:length(plotGroups{ind})-1
                    currStartTime = segmentInfo.timeEnd(plotGroups{ind}(ind_seg));
                    
                    switch ind_seg
                        case length(plotGroups{ind})
                            
                            switch ind
                                case length(plotGroups)
                                    % last entry in the array entirely
                                    currEndTime = feature_full.t(end);
                                    
                                otherwise
                                    % first entry in the next group
                                    currEndTime = segmentInfo.timeStart(plotGroups{ind+1}(1));
                            end
                            
                        otherwise
                            % next entry in this group
                            currEndTime = segmentInfo.timeStart(plotGroups{ind}(ind_seg+1));
                    end
                    
                    
                    [startVal, startInd] = findClosestValue(currStartTime, feature_full.t);
                    [endVal, endInd] = findClosestValue(currEndTime, feature_full.t);
                    
                    indsToUse = [indsToUse startInd:endInd];
                 end
        end
        
        meanWeight_group(ind, :) = mean(avgWeightArray_belowThres(indsToUse, :), 1);
        stdWeight_group(ind, :) = std(avgWeightArray_belowThres(indsToUse, :), 1);
        
        meanResnrom_group(ind, :) = mean(resnromAll_lsqlin_const_minRMSE_array_belowThres(indsToUse, :), 1);
        stdResnrom_group(ind, :) = std(resnromAll_lsqlin_const_minRMSE_array_belowThres(indsToUse, :), 1);
    else
        meanWeight_group(ind, :) = zeros(1, size(avgWeightArray_belowThres, 2));
        stdWeight_group(ind, :) = zeros(1, size(avgWeightArray_belowThres, 2));
        
        meanResnrom_group(ind, :) = zeros(1, size(resnromAll_lsqlin_const_minRMSE_array_belowThres, 2));
        stdResnrom_group(ind, :) = zeros(1, size(resnromAll_lsqlin_const_minRMSE_array_belowThres, 2));
    end
end

figure(h);

% pull out the top values
maxAvgWeight = mean(avgWeightArray_belowThres);
[~, maxInd1] = max(maxAvgWeight);
maxAvgWeight(maxInd1) = 0;
[~, maxInd2] = max(maxAvgWeight);
maxAvgWeight(maxInd2) = 0;
[~, maxInd3] = max(maxAvgWeight);
maxAvgWeight(maxInd3) = 0;

top3ind = sort([maxInd1, maxInd2, maxInd3]);

ax(1) = subplot(2, 1, 1); hold on
% title('Averaged Normalized Cost Weights');

x = 1:length(plotGroups);
if length(x) == 1
    barwitherr(stdWeight_group(:, top3ind), meanWeight_group(:, top3ind));
else
    barwitherr(stdWeight_group(:, top3ind), 1:length(plotGroups), meanWeight_group(:, top3ind));
end
legend(cost_function_names(top3ind));
title('Cost function');

% bar(feature_full.t, avgWeightArray_belowThres, 'stacked'); 
% shading flat
% legend(cost_function_names);
% % xlabel('Time [s]');
% ylabel('Normalized Cost Weights');
% 
% ylim([-0.05 1.05])

ax(2) = subplot(2, 1, 2); hold on

x = 1:length(plotGroups);
if length(x) == 1
    barwitherr(stdResnrom_group(:, :), meanResnrom_group(:, :));
else
    barwitherr(stdResnrom_group(:, :), 1:length(plotGroups), meanResnrom_group(:, :));
end

title('Resnorm');

% % title('Averaged Residual Norm');
% plot(xlim, [constThres constThres], 'k');
% xlabel('Time [s]');
% ylabel('Averaged Residual Norm');
% 
% % upperbound = min([constThres*5 max(resnromAll_lsqlin_const_minRMSE_array_belowThres)]);
% 
% ylim([0 constThres*5]);
% 
linkaxes(ax, 'x');
% xlim([t_rmse(1) - 0.5 t_rmse(end) + 0.5]);

% also generate a correlation plot
for ind_costf = 1:length(cost_function_names)
    x_val = [1:length(plotGroups)]';
    y1_val = meanWeight_group(:, ind_costf);
    p1 = polyfit(x_val, y1_val, 1);
    f1 = polyval(p1, x_val);
    [rmse_b1, rmse_y1, rmse_Rsq1] = calcRegression(x_val, y1_val);
     
    y2_val = stdWeight_group(:, ind_costf);
    p2 = polyfit(x_val, y2_val, 1);
    f2 = polyval(p2, x_val);
    [rmse_b2, rmse_y2, rmse_Rsq2] = calcRegression(x_val, y2_val);
    
%     h28 = figure;
%     subplot(211);
%     plot(x_val, y1_val, 'o'); hold on
%     plot(x_val, f1, 'b-');
%     title([cost_function_names{ind_costf} ' weight mean R^2 = ' num2str(rmse_Rsq1)]);
%     
%     subplot(212);
%     plot(x_val, y2_val, 'o'); hold on
%     plot(x_val, f2, 'b-');
%     title([cost_function_names{ind_costf} ' weight stddev R^2 = ' num2str(rmse_Rsq2)]);
%     
%     saveas(h28, fullfile(outputPath, [currInstName  num2str(ind_costf) '_corr_'  '.fig']));
%     saveas(h28, fullfile(outputPath, [currInstName num2str(ind_costf) '_corr_'  '.png']));
%     
%     close(h28);
    
    rsq_mean(ind_costf) = rmse_Rsq1;
    rsq_std(ind_costf) = rmse_Rsq2;
end

    x_val = [1:length(plotGroups)]';
    y1_val = meanResnrom_group(:, :);
    p1 = polyfit(x_val, y1_val, 1);
    f1 = polyval(p1, x_val);
    [rmse_b1, rmse_y1, rsq_mean(end+1)] = calcRegression(x_val, y1_val);
     
    y2_val = stdResnrom_group(:, :);
    p2 = polyfit(x_val, y2_val, 1);
    f2 = polyval(p2, x_val);
    [rmse_b2, rmse_y2, rsq_std(end+1)] = calcRegression(x_val, y2_val);
    
%       h28 = figure;
%     subplot(211);
%     plot(x_val, y1_val, 'o'); hold on
%     plot(x_val, f1, 'b-');
%     title([ 'resnorm weight mean R^2 = ' num2str(rmse_Rsq1)]);
%     
%     subplot(212);
%     plot(x_val, y2_val, 'o'); hold on
%     plot(x_val, f2, 'b-');
%     title([ 'resnorm weight stddev R^2 = ' num2str(rmse_Rsq2)]);
%     
%     saveas(h28, fullfile(outputPath, [currInstName  num2str(ind_costf+1) '_corr_'  '.fig']));
%     saveas(h28, fullfile(outputPath, [currInstName num2str(ind_costf+1) '_corr_'  '.png']));
%     
%     close(h28);