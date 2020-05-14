function plotCalc_applyThreshold(output_inverse, indToUse_window, cost_function_names, constThres, noWinMode, plotType)
    % constThres - the threshold to reject high resnorms
    % noWinMode - what to do if a given window has no windows attached 
    
    
% constThres = 10e-3;    % check1threshold 0.1e-3 1e-3 (1 DOF), 10e-3 (3
% DOF), for offset mode
diffThreshold = 0;  % check2threshold 1e-3
negConstThreshold = -0.05;

% individual cost functions
windowCount = size(indToUse_window, 1);
ccostCount = length(cost_function_names);

cAll_const_plot = [];
cAll_unconst_plot = [];

resnormAll_sliding_const_plot = [];
resnormAll_sliding_unconst_plot = [];

q_blended = [];
windowed_rmse_belowThreshold = [];

for ind_windowCount = 1:windowCount
    % pull out the offset version of each data
    switch param.segment_only_windows
        case 'none'
            frontEdgeInd = indToUse_window(ind_windowCount, 2)-param.win_shift+1;
            if frontEdgeInd < 1
                frontEdgeInd = 1;
            end
            
        otherwise
            frontEdgeInd = indToUse_window(ind_windowCount, 1);
    end
    
    ind_toUseT = frontEdgeInd:indToUse_window(ind_windowCount, 2); % the time length covered is the same as the window inserted
    lenToUseT = length(ind_toUseT);
    
    switch plotType
        case 'weight'
            c_recovered_curr = output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).c_recovered;
            
        case 'ratio'
            c_recovered_curr = output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).J_out_contrib;
    end
    
    c_recovered_curr = c_recovered_curr / sum(c_recovered_curr);
    
    c_unconst_curr = output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).c_out_pivot_unconst;
    resnormAll_sliding_const_plot  =  [resnormAll_sliding_const_plot;    repmat([output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).resnorm_lsqlin],         [lenToUseT 1])];
    resnormAll_sliding_unconst_plot = [resnormAll_sliding_unconst_plot;  repmat([output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).resnorm_lsqlin_unconst], [lenToUseT 1])];
    cAll_const_plot =     [cAll_const_plot;   repmat(c_recovered_curr,             [lenToUseT 1])];
    cAll_unconst_plot =   [cAll_unconst_plot; repmat(c_unconst_curr,               [lenToUseT 1])];
    
    resnorm_curr = output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).resnorm_lsqlin;
    rmse_temp = calc_rmse(rmse_fct, feature_win_save{ind_windowCount}, feature_recon_local{ind_windowCount}{minRmseIndArray(ind_windowCount)}, param);
    windowed_rmse_all(ind_windowCount) = [rmse_temp.q];
    
    if resnorm_curr <= constThres
        % now calculate windowed RMSE
        t_belowThreshold{ind_windowCount} = t_recon_plot_array{ind_windowCount};
        q_belowThreshold{ind_windowCount} = q_recon_plot_array{ind_windowCount};
        feat_belowThreshold{ind_windowCount} = feature_recon_local{ind_windowCount}{minRmseIndArray(ind_windowCount)};
        
        windowed_rmse_belowThreshold = [windowed_rmse_belowThreshold rmse_temp.q];
    else
        t_belowThreshold{ind_windowCount} = {};
        q_belowThreshold{ind_windowCount} = {};
        feat_belowThreshold{ind_windowCount} = {};
    end
end 

t_belowThreshold_combined_temp = horzcat(t_belowThreshold{:}); %horzcat(J_inverse{:});
q_belowThreshold_combined_temp = horzcat(q_belowThreshold{:});

if iscell(t_belowThreshold_combined_temp)
    t_belowThreshold_combined = horzcat(t_belowThreshold_combined_temp{:});
    q_belowThreshold_combined = horzcat(q_belowThreshold_combined_temp{:});
else
    t_belowThreshold_combined = t_belowThreshold_combined_temp;
    q_belowThreshold_combined = q_belowThreshold_combined_temp;
end

rmsEntryToUse = zeros(1, length(feature_full.t));

% now for the average version. at each timestep, which window corresponds?
for ind_tCount = 1:length(feature_full.t)  
    % find all the corresponding timesteps
    startCheck = ind_tCount >= indToUse_window(:, 1);
    endCheck =                 indToUse_window(:, 2) >= ind_tCount;
    
    % did it appear in both lists? 
    bothCheck = and(startCheck, endCheck);
    valToUse = find(bothCheck);
    
    if ~isempty(valToUse)
        c_recovered_curr = [];
        c_negConst_recovered_curr = [];
        resnormAll_lsqlin_const_RMSE_curr =   [];
        resnormAll_lsqlin_unconst_RMSE_curr = [];
        
        for ind_check = 1:length(valToUse)
            curr_valToUse = valToUse(ind_check);
            switch plotType
                case 'weight'
                    c_recovered_curr = [c_recovered_curr; output_inverse{curr_valToUse}(minRmseIndArray(curr_valToUse)).c_recovered];
                        
                case 'ratio'
                    c_recovered_curr = [c_recovered_curr; output_inverse{curr_valToUse}(minRmseIndArray(curr_valToUse)).J_out_contrib];
            end
            
            c_negConst_recovered_curr = [c_negConst_recovered_curr; output_inverse{curr_valToUse}(minRmseIndArray(curr_valToUse)).c_out_pivot_unconst];
         
            resnormAll_lsqlin_const_RMSE_curr =   [resnormAll_lsqlin_const_RMSE_curr;   output_inverse{curr_valToUse}(minRmseIndArray(curr_valToUse)).resnorm_lsqlin];
            resnormAll_lsqlin_unconst_RMSE_curr = [resnormAll_lsqlin_unconst_RMSE_curr; output_inverse{curr_valToUse}(minRmseIndArray(curr_valToUse)).resnorm_lsqlin_unconst];
        end
        
        c_temp1 = mean(c_recovered_curr, 1);
        c_temp2 = mean(c_negConst_recovered_curr, 1);
        res_temp1 = mean(resnormAll_lsqlin_const_RMSE_curr, 1);
        res_temp2 = mean(resnormAll_lsqlin_unconst_RMSE_curr, 1);
        
        c_temp1 = c_temp1/sum(c_temp1);
        c_temp2 = c_temp2/sum(c_temp2);
%         res_temp1 = res_temp1/sum(res_temp1);
%         res_temp2 = res_temp2/sum(res_temp2);
        
        avgWeightArray(ind_tCount, :) = c_temp1;
        avgWeightNegConstArray(ind_tCount, :) = c_temp2;
        resnromAll_lsqlin_const_minRMSE_array(ind_tCount, :) = res_temp1;
        resnromAll_lsqlin_unconst_minRMSE_array(ind_tCount, :) = res_temp2;
        
        % remove the entries that are above the threshold
        keepInd = find(resnormAll_lsqlin_const_RMSE_curr < constThres);
        
        if isempty(keepInd)
            switch noWinMode
                case 'blank'
                    c_recover_belowThres = zeros(1, ccostCount);
                    resnorm_belowThres = 0;
                case 'previous'
                    if ind_tCount == 1
                        c_recover_belowThres = zeros(1, ccostCount);
                        resnorm_belowThres = 0;
                    else
                        c_recover_belowThres = avgWeightArray_belowThres(ind_tCount - 1, :);
                        resnorm_belowThres = resnromAll_lsqlin_const_minRMSE_array_belowThres(ind_tCount - 1, :);
                    end
            end
        else
            c_temp = c_recovered_curr(keepInd, :);
            res_temp = resnormAll_lsqlin_const_RMSE_curr(keepInd, :);
            
            c_temp = mean(c_temp, 1);
            res_temp = mean(res_temp, 1);
            
            c_recover_belowThres = c_temp / (sum(c_temp));
            resnorm_belowThres = res_temp;
        end
        
        avgWeightArray_belowThres(ind_tCount, :) = c_recover_belowThres;
        resnromAll_lsqlin_const_minRMSE_array_belowThres(ind_tCount, :) = resnorm_belowThres;
        resnromAll_lsqlin_const_minRMSE_array_belowThresPass(ind_tCount) = length(keepInd);
    else
        avgWeightArray(ind_tCount, :) = zeros(1, ccostCount);
        avgWeightNegConstArray(ind_tCount, :) = zeros(1, ccostCount);
        
        resnromAll_lsqlin_const_minRMSE_array(ind_tCount, :) = 0;
        resnromAll_lsqlin_unconst_minRMSE_array(ind_tCount, :) = 0;
        
        avgWeightArray_belowThres(ind_tCount, :) = zeros(1, ccostCount);
        resnromAll_lsqlin_const_minRMSE_array_belowThres(ind_tCount, :) = 0;
    end
    
     % now blend together the q over all the pass windows
     currT = feature_full.t(ind_tCount);     
     findTimeStep = find(t_belowThreshold_combined == currT);
     if ~isempty(findTimeStep) 
         % pull in existing values
         findQ = q_belowThreshold_combined(:, findTimeStep);
         q_blended(:, ind_tCount) = mean(findQ, 2);
         rmsEntryToUse(ind_tCount) = 1;
     elseif isempty(findTimeStep) && isempty(q_blended)
         % no value defined for that entry, and nothing to set as previous
         q_blended(:, ind_tCount) = zeros(size(q_belowThreshold_combined, 1), 1);
         rmsEntryToUse(ind_tCount) = 0;
     else
         switch noWinMode
             case 'blank'
                 q_blended(:, ind_tCount) = zeros(size(q_belowThreshold_combined, 1), 1);
                 rmsEntryToUse(ind_tCount) = 0;
                 
             case 'previous'
                 q_blended(:, ind_tCount) = q_blended(:, ind_tCount-1);
                 rmsEntryToUse(ind_tCount) = 1;
         end
     end
end

% indToUseStart = indToUse_window(1, 1); 
% indToUseEnd = indToUse_window(end, 2);
blended_rmse_belowThreshold = rmse_fct(feature_full.q(:, find(rmsEntryToUse)), q_blended(:, find(rmsEntryToUse)), []);

ccost_full_plot = ccost_full ./ repmat(sum(ccost_full, 2), 1, size(ccost_full, 2));

% resnorm rejection criteria
meanconstResnrom = mean(resnormAll_sliding_const_plot);
meanunconstResnorm = mean(resnormAll_sliding_unconst_plot);

meanconstResnrom_averageWindow  = mean(resnromAll_lsqlin_const_minRMSE_array);
meanunconstResnorm_averageWindow  = mean(resnormAll_sliding_unconst_plot);

pts_aboveConstThres = resnormAll_sliding_const_plot >= constThres;
pts_underConstThres = resnormAll_sliding_const_plot <  constThres;

pts_aboveDiffThres = resnormAll_sliding_const_plot - resnormAll_sliding_unconst_plot >= diffThreshold;
pts_underDiffThres = resnormAll_sliding_const_plot - resnormAll_sliding_unconst_plot <  diffThreshold;

pts_allow = or(pts_underConstThres, pts_underDiffThres);
pts_deny = setxor(1:length(pts_allow), find(pts_allow));
% pts_deny =  pts_aboveConstThres && pts_aboveDiffThres;

pts_aboveConstThres_averageWindow = resnromAll_lsqlin_const_minRMSE_array >= constThres;
pts_underConstThres_averageWindow = resnromAll_lsqlin_const_minRMSE_array <  constThres;

minWeightsUnconst = min(cAll_unconst_plot, [], 2);
lowestValInUnconst = minWeightsUnconst > negConstThreshold;

% now we can do our metric counting
resnormPass1 = sum(pts_underConstThres);
resnormPass2 = sum(pts_underDiffThres);
resnormPass3 = sum(lowestValInUnconst);
resnormPass4 = sum(pts_underConstThres_averageWindow);

weightsActive = sum(avgWeightArray_belowThres);
weightsMean = mean(avgWeightArray_belowThres);

% set them up in structs so we can tally them 
rmse_report.resnormPass1 = resnormPass1;
rmse_report.resnormPass2 = resnormPass2;
rmse_report.resnormPass3 = resnormPass3;
rmse_report.resnormPass4 = resnormPass4;

rmse_report.weightsActive = weightsActive;
rmse_report.weightsTotal = length(feature_full.t) ;
rmse_report.weightsMean = weightsMean;

rmse_report.blended_rmse_belowThreshold = blended_rmse_belowThreshold;
rmse_report.windowed_rmse_belowThreshold_mean = mean(windowed_rmse_belowThreshold);
rmse_report.windowed_rmse_belowThreshold_std = std(windowed_rmse_belowThreshold);
end