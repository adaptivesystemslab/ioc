% plot data
[~, uniqueInd] = unique(t_rmse);
ax = [];

titles = {'full window in static', 'full window in motion', 'full window in both'};

% calculate the mean and range of the contribution ratios for all rest
% (not inside a segment window) and segment (all motions within a segment
% window) 
ratioAll_end{1} = []; % the window ends stationary
ratioAll_end{2} = []; % the window ends in a motion
ratioAll_all{1} = []; % the whole window is not inside a motion
ratioAll_all{2} = []; % the whole window is inside a motion
ratioAll_all{3} = []; % everything else

for ind_windowCount = 1:windowCount
    t_windowCount = feature_full.t(:, indToUse_window(ind_windowCount, 1):indToUse_window(ind_windowCount, 2));

    [motionEnd, motionWhole, ratioAll_end, ratioAll_all] = ...
        checkStandingVSMovement(feature_full, indToUse_window(ind_windowCount, 1):indToUse_window(ind_windowCount, 2), ...
        segmentInfo, ratioRMSE{ind_windowCount}, ratioAll_end, ratioAll_all);

end

h4 = figure;
ratioMean = [];
ratioSD = [];
for ii = 1:3
%     ax(ii) = subplot(3, 1, ii); hold on
    
    % calculate the mean and std
    if isempty(ratioAll_all{ii})
        ratioMean(ii, :) = zeros(1, length(cost_function_names));
        ratioSD(ii, :) = zeros(1, length(cost_function_names));
    else
        ratioMean(ii, :) = mean(ratioAll_all{ii}, 1);
        ratioSD(ii, :) = std(ratioAll_all{ii}, 0, 1);    
    end
end

h = barwitherr(ratioSD', ratioMean');% Plot with errorbars
legend(titles);
% xlim([0.5 length(cost_function_names)+0.5]);
set(gca,'XTickLabel',cost_function_names);
title('Contrib ratio in stationary/movement');
