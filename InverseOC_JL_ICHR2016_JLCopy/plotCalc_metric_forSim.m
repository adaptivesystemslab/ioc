plotType = 'weight';

constThres = returnThreshold(param);

% constThres = 10e-3;    % check1threshold 0.1e-3 1e-3 (1 DOF), 10e-3 (3
% DOF), for offset mode
diffThreshold = 0;  % check2threshold 1e-3
negConstThreshold = -0.05;

% individual cost functions
windowCount = size(indToUse_window, 1);
ccostCount = length(cost_function_names);

cAll_plot = [];
cAllUnconst_plot = [];

resnormAll_lsqlin_const_RMSE_plot = [];
resnormAll_lsqlin_unconst_RMSE_plot = [];

for ind_windowCount = 1:windowCount
    % pull out the offset version of each data
    switch plotType
        case 'weight'
            c_recovered_curr = output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).c_recovered;
                 
        case 'ratio'
            c_recovered_curr = output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).J_out_contrib;
    end
    
    c_unconst_curr = output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).c_out_pivot_unconst;
    resnormAll_lsqlin_const_RMSE_plot =  [resnormAll_lsqlin_const_RMSE_plot;   repmat([output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).resnorm_lsqlin], [lenToUseT 1])];
    resnormAll_lsqlin_unconst_RMSE_plot = [resnormAll_lsqlin_unconst_RMSE_plot; repmat([output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).resnorm_lsqlin_unconst], [lenToUseT 1])];
    
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
    
    cAll_plot =   [cAll_plot;   repmat(c_recovered_curr,                           [lenToUseT 1])];
    cAllUnconst_plot =   [cAllUnconst_plot;   repmat(c_unconst_curr,                           [lenToUseT 1])];
end

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
            
        avgWeightArray(ind_tCount, :) = mean(c_recovered_curr, 1);
        avgWeightNegConstArray(ind_tCount, :) = mean(c_negConst_recovered_curr, 1);
        
        resnromAll_lsqlin_const_minRMSE_array(ind_tCount, :) = mean(resnormAll_lsqlin_const_RMSE_curr, 1);
        resnromAll_lsqlin_unconst_minRMSE_array(ind_tCount, :) = mean(resnormAll_lsqlin_unconst_RMSE_curr, 1);
        
        % remove the entries that are above the threshold
        keepInd = find(resnormAll_lsqlin_const_RMSE_curr < constThres);
        
           
        if isempty(keepInd)
            c_recover_belowThres = zeros(1, ccostCount);
            resnorm_belowThres = 0;
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
        
        if length(resnromAll_lsqlin_const_minRMSE_array_belowThres) ~= ind_tCount
            slkfaj = 1; 
        end
        
    else
        avgWeightArray(ind_tCount, :) = zeros(1, ccostCount);
        avgWeightNegConstArray(ind_tCount, :) = zeros(1, ccostCount);
        
        resnromAll_lsqlin_const_minRMSE_array(ind_tCount, :) = 0;
        resnromAll_lsqlin_unconst_minRMSE_array(ind_tCount, :) = 0;
        
        avgWeightArray_belowThres(ind_tCount, :) = zeros(1, ccostCount);
        resnromAll_lsqlin_const_minRMSE_array_belowThres(ind_tCount, :) = 0;
    end
end


% resnorm rejection criteria
meanconstResnrom = mean(resnormAll_lsqlin_const_RMSE_plot);
meanunconstResnorm = mean(resnormAll_lsqlin_unconst_RMSE_plot);

meanconstResnrom_averageWindow  = mean(resnromAll_lsqlin_const_minRMSE_array);
meanunconstResnorm_averageWindow  = mean(resnormAll_lsqlin_unconst_RMSE_plot);

pts_aboveConstThres = resnormAll_lsqlin_const_RMSE_plot >= constThres;
pts_underConstThres = resnormAll_lsqlin_const_RMSE_plot <  constThres;

pts_aboveDiffThres = resnormAll_lsqlin_const_RMSE_plot - resnormAll_lsqlin_unconst_RMSE_plot >= diffThreshold;
pts_underDiffThres = resnormAll_lsqlin_const_RMSE_plot - resnormAll_lsqlin_unconst_RMSE_plot <  diffThreshold;

pts_allow = or(pts_underConstThres, pts_underDiffThres);
pts_deny = setxor(1:length(pts_allow), find(pts_allow));
% pts_deny =  pts_aboveConstThres && pts_aboveDiffThres;

pts_aboveConstThres_averageWindow = resnromAll_lsqlin_const_minRMSE_array >= constThres;
pts_underConstThres_averageWindow = resnromAll_lsqlin_const_minRMSE_array <  constThres;

minWeightsUnconst = min(cAllUnconst_plot, [], 2);
lowestValInUnconst = minWeightsUnconst > negConstThreshold;

% now we can do our metric counting
resnormPass1 = sum(pts_underConstThres);
resnormPass2 = sum(pts_underDiffThres);
resnormPass3 = sum(lowestValInUnconst);
resnormPass4 = sum(pts_underConstThres_averageWindow);

weightsActive = sum(avgWeightArray > 0);
weightsMean = mean(avgWeightArray);

% set them up in structs so we can tally them 
rmse_report.resnormPass1 = resnormPass1;
rmse_report.resnormPass2 = resnormPass2;
rmse_report.resnormPass3 = resnormPass3;
rmse_report.resnormPass4 = resnormPass4;

rmse_report.weightsActive = weightsActive;
rmse_report.weightsMean = weightsMean;