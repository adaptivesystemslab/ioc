% constThres = returnThreshold(param);

% constThresMultiplier = 1;
constThresMultiplier = thresholdMultiplier;

switch param.dataset
%     case {'healthy1ankle'}
% %         constThres = 2.21*0.95*constThresMultiplier; % squat ankle
%         constThres = 3.8*0.95*constThresMultiplier; % squat ankle
%         
%     case {'healthy1hip', 'sim'}
% %         constThres = 10.8*0.95*constThresMultiplier; % squat hip
%         constThres = 5.75*0.95*constThresMultiplier; % squat hip
%         
%     case {'squats_tuat_2011'}
% %         constThres = 0.0923*0.95*constThresMultiplier;
%         constThres = 17.4*0.95*constThresMultiplier;
%         
%     case {'squats_tuat_2015'}
% %         constThres = 0.736*0.95*constThresMultiplier; % squat ufri
%         constThres = 2.32*0.95*constThresMultiplier; % squat ufri
        
    otherwise
        constThres = 1e9;
end
    
% constThres = 0.24; % squat ufri
% constThres = 0.73; % squat ankle

% constThres = 1e9;
% svd_ratio_threshold = 0.99;

specialDof = [1:length(cost_function_names)];

% no window mode? 
noWinMode = 'blank';

% constThres = 10e-3;    % check1threshold 0.1e-3 1e-3 (1 DOF), 10e-3 (3
% DOF), for offset mode
diffThreshold = 0;  % check2threshold 1e-3
negConstThreshold = -0.05;

indName = 0;
for i = 1:126
    indName = indName + 1;
    lambda_name{indName} = ['lambda' num2str(i)];
end

% individual cost functions
windowCount = size(indToUse_window, 1);
ccostCount = length(cost_function_names);
cconstCount = length(lambda_name);

t_rmse = [];
cAll_plot = [];
% cAllUnconst_plot = [];

resnormAll_lsqlin_const_RMSE_plot = [];
% resnormAll_lsqlin_unconst_RMSE_plot = [];
% correlationArray_special{ccostCount} = [];
recoverWeightAll{ccostCount} = [];
resnorm_special{ccostCount} = [];
rmse_special{ccostCount} = [];
rankGrad_special{ccostCount+1} = [];
svd_special{ccostCount+1} = [];
rankA_special{ccostCount} = [];
dofsUsedInIOC_plot = [];

resnorm_min = [];
rmse_min = [];

allDof_plot = [];
minDof_plot = [];

q_blended = [];
dq_blended = [];
ddq_blended = [];
windowed_rmse_belowThreshold = [];
windowed_rmse_belowThreshold_dq = [];
windowed_rmse_belowThreshold_ddq = [];

 lambda_array = [];

for ind_windowCount = 1:windowCount
    % pull out the offset version of each data
    minDof = minRmseIndArray(ind_windowCount);
    
     switch param.segment_only_windows
        case 'none'
            frontEdgeInd = indToUse_window(ind_windowCount, 2)-param.win_shift+1;
            if frontEdgeInd < 1
                frontEdgeInd = 1;
            end
            
        otherwise
            frontEdgeInd = indToUse_window(ind_windowCount, 1);
     end
    
     endEdgeInd = indToUse_window(ind_windowCount, 2);
     
%      if endEdgeInd > length(feature_full.t)
%          endEdgeInd = length(feature_full.t);
%      end
    
    ind_toUseT = frontEdgeInd:endEdgeInd; % the time length covered is the same as the window inserted
    lenToUseT = length(ind_toUseT);
    
    switch plotType
        case 'weight'
            c_recovered_curr = output_inverse{ind_windowCount}(minDof).c_recovered;
                 
        case 'ratio'
            c_recovered_curr = output_inverse{ind_windowCount}(minDof).J_out_contrib;
    end
    
    c_recovered_curr = c_recovered_curr / sum(c_recovered_curr);
    
%     c_unconst_curr = output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).c_out_pivot_unconst;
    resnormAll_lsqlin_const_RMSE_plot =  [resnormAll_lsqlin_const_RMSE_plot;   repmat([output_inverse{ind_windowCount}(minDof).resnorm_lsqlin], [lenToUseT 1])];
%     resnormAll_lsqlin_unconst_RMSE_plot = [resnormAll_lsqlin_unconst_RMSE_plot; repmat([output_inverse{ind_windowCount}(minDof).resnorm_lsqlin_unconst], [lenToUseT 1])];
        
    for j = 1:length(specialDof)
        localSpecialDof = specialDof(j);
%         currSpecial = output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).corrValMatrix_All(:, localSpecialDof);
%         currCorrSize = size(currSpecial);
%         correlationArray_special{j} = [correlationArray_special{j}; repmat(reshape(currSpecial, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];

        currSpecial = output_inverse{ind_windowCount}(localSpecialDof).c_recovered;
        currSpecial = currSpecial / sum(currSpecial);
        currCorrSize = size(currSpecial);
        recoverWeightAll{j} = [recoverWeightAll{j}; repmat(reshape(currSpecial, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];

        currSpecial = output_inverse{ind_windowCount}(localSpecialDof).resnorm_lsqlin;
        currCorrSize = size(currSpecial);
        resnorm_special{j} = [resnorm_special{j}; repmat(reshape(currSpecial, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];

        currSpecial = output_inverse{ind_windowCount}(localSpecialDof).rmse;
        currCorrSize = size(currSpecial);
        rmse_special{j} = [rmse_special{j}; repmat(reshape(currSpecial, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];

% %          currSpecial = output_inverse{ind_windowCount}(localSpecialDof).eig_cutoff;
%         currSpecial = output_inverse{ind_windowCount}(localSpecialDof).J_coeff_individual_out(:, [1:(localSpecialDof-1) (localSpecialDof+1):end]);
% %         currSpecial = output_inverse{ind_windowCount}(localSpecialDof).A;
%         currSpecial = rank(currSpecial);
%         currCorrSize = size(currSpecial);
%         rankGrad_special{j} = [rankGrad_special{j}; repmat(reshape(currSpecial, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];
        
        currSpecial = output_inverse{ind_windowCount}(localSpecialDof).eig_cutoff;
        svd_special{j} = [svd_special{j}; repmat(reshape(currSpecial, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];
        
% %         currSpecial = output_inverse{ind_windowCount}(localSpecialDof).J_coeff_individual_out(:, [1:(localSpecialDof-1) (localSpecialDof+1):end]);
%         currSpecial = output_inverse{ind_windowCount}(localSpecialDof).A;
%         currSpecial = rank(currSpecial);
%         currCorrSize = size(currSpecial);
%         rankA_special{j} = [rankA_special{j}; repmat(reshape(currSpecial, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];
    end
    
%     currSpecial = output_inverse{ind_windowCount}(minDof).J_coeff_individual_out;
%     currSpecial = rank(currSpecial);
%     currCorrSize = size(currSpecial);
%     rankGrad_special{j+1} = [rankGrad_special{j+1}; repmat(reshape(currSpecial, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];
        
    svd_special{j+1} = [svd_special{j+1}; repmat(reshape(currSpecial, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];
       
    currSpecial = output_inverse{ind_windowCount}(minDof).rref_out_full;
    rref_rank(ind_windowCount, 1) = output_inverse{ind_windowCount}(minDof).eig_cutoff_full;
    rref_rank(ind_windowCount, 2) = feature_full.t(:, ind_toUseT(1));
    rref_out{ind_windowCount} = currSpecial;
    
    currSpecial = output_inverse{ind_windowCount}(minDof).resnorm_lsqlin;
    currCorrSize = size(currSpecial);
    resnorm_min = [resnorm_min; repmat(reshape(currSpecial, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];

    currSpecial = output_inverse{ind_windowCount}(minDof).rmse;
    currCorrSize = size(currSpecial);
    rmse_min = [rmse_min; repmat(reshape(currSpecial, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];
    
    currSpecial = output_inverse{ind_windowCount}(minDof).cost_functions_used_correlation;
    zeroMat = zeros(size(cost_function_names));
    zeroMat(currSpecial) = currSpecial;
    currCorrSize = size(zeroMat);
    allDof_plot = [allDof_plot; repmat(reshape(zeroMat, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];    
   
    currSpecial = output_inverse{ind_windowCount}(minDof).lambda';
    currCorrSize = size(currSpecial);
    lambda_array = [lambda_array; repmat(reshape(currSpecial, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];
    
    currSpecial = minDof;
    currCorrSize = size(currSpecial);
    minDof_plot = [minDof_plot; repmat(reshape(currSpecial, 1, currCorrSize(1)*currCorrSize(2)), [lenToUseT 1])];
    
    t_rmse = [t_rmse feature_full.t(:, ind_toUseT)]; % cover the increment between the previous window and this one
    cAll_plot =   [cAll_plot;   repmat(c_recovered_curr,                           [lenToUseT 1])];
%     cAllUnconst_plot =   [cAllUnconst_plot;   repmat(c_unconst_curr,               [lenToUseT 1])];
    
    currResNorm = output_inverse{ind_windowCount}(minRmseIndArray(ind_windowCount)).resnorm_lsqlin;
    rmse_temp = calc_rmse(rmse_fct, feature_win_save{ind_windowCount}, feature_recon_local{ind_windowCount}{minRmseIndArray(ind_windowCount)}, param);
    windowed_rmse_all(ind_windowCount) = [rmse_temp.q];
    
    if currResNorm <= constThres
        % now calculate windowed RMSE
        t_belowThreshold{ind_windowCount} = t_recon_plot_array{ind_windowCount};
        q_belowThreshold{ind_windowCount} = q_recon_plot_array{ind_windowCount};
        dq_belowThreshold{ind_windowCount} = dq_recon_plot_array{ind_windowCount};
        ddq_belowThreshold{ind_windowCount} = ddq_recon_plot_array{ind_windowCount};
        feat_belowThreshold{ind_windowCount} = feature_recon_local{ind_windowCount}{minRmseIndArray(ind_windowCount)};
        
        windowed_rmse_belowThreshold = [windowed_rmse_belowThreshold rmse_temp.q];
        windowed_rmse_belowThreshold_dq = [windowed_rmse_belowThreshold_dq rmse_temp.dq];
        windowed_rmse_belowThreshold_ddq = [windowed_rmse_belowThreshold_ddq rmse_temp.ddq];
    else
        t_belowThreshold{ind_windowCount} = {};
        q_belowThreshold{ind_windowCount} = {};
        dq_belowThreshold{ind_windowCount} = {};
        ddq_belowThreshold{ind_windowCount} = {};
        feat_belowThreshold{ind_windowCount} = {};
    end
end 

% lambda_array_belowThreshold = lambda_array;

% figure;
% plot(correlationArray_special(:, 1:4));
% figure;
% plot(correlationArray_special(:, 5:8));
% figure;
% plot(correlationArray_special(:, 9:end));


t_belowThreshold_combined_temp = horzcat(t_belowThreshold{:}); %horzcat(J_inverse{:});
q_belowThreshold_combined_temp = horzcat(q_belowThreshold{:});
dq_belowThreshold_combined_temp = horzcat(dq_belowThreshold{:});
ddq_belowThreshold_combined_temp = horzcat(ddq_belowThreshold{:});

if iscell(t_belowThreshold_combined_temp)
    t_belowThreshold_combined = horzcat(t_belowThreshold_combined_temp{:});
    q_belowThreshold_combined = horzcat(q_belowThreshold_combined_temp{:});
    dq_belowThreshold_combined = horzcat(dq_belowThreshold_combined_temp{:});
    ddq_belowThreshold_combined = horzcat(ddq_belowThreshold_combined_temp{:});
else
    t_belowThreshold_combined = t_belowThreshold_combined_temp;
    q_belowThreshold_combined = q_belowThreshold_combined_temp;
    dq_belowThreshold_combined = dq_belowThreshold_combined_temp;
    ddq_belowThreshold_combined = ddq_belowThreshold_combined_temp;
end

rmsEntryToUse = zeros(1, length(feature_full.t));
[withinSegmentCheck] = checkWithinSegment(feature_full, feature_full.t, segmentInfo);

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
        J_recovered_curr = [];
        J_recovered_abs_curr = [];
        lambda_recovered_curr = [];
%         c_negConst_recovered_curr = [];
        resnormAll_lsqlin_const_RMSE_curr =   [];
        rmseAll_lsqlin_const_RMSE_curr =   [];
%         resnormAll_lsqlin_unconst_RMSE_curr = [];
        
        for ind_check = 1:length(valToUse)
            curr_valToUse = valToUse(ind_check);
            
            c_recovered_curr = [c_recovered_curr; output_inverse{curr_valToUse}(minRmseIndArray(curr_valToUse)).c_recovered];
            J_recovered_curr = [J_recovered_curr; output_inverse{curr_valToUse}(minRmseIndArray(curr_valToUse)).J_out_contrib];
            J_recovered_abs_curr = [J_recovered_abs_curr; output_inverse{curr_valToUse}(minRmseIndArray(curr_valToUse)).J_out_array];
            lambda_recovered_curr = [lambda_recovered_curr; output_inverse{curr_valToUse}(minRmseIndArray(curr_valToUse)).lambda'];
           
%             c_negConst_recovered_curr = [c_negConst_recovered_curr; output_inverse{curr_valToUse}(minRmseIndArray(curr_valToUse)).c_out_pivot_unconst];
         
            resnormAll_lsqlin_const_RMSE_curr =   [resnormAll_lsqlin_const_RMSE_curr;   output_inverse{curr_valToUse}(minRmseIndArray(curr_valToUse)).resnorm_lsqlin];
            rmseAll_lsqlin_const_RMSE_curr =   [rmseAll_lsqlin_const_RMSE_curr;   output_inverse{curr_valToUse}(minRmseIndArray(curr_valToUse)).rmse];
%             resnormAll_lsqlin_unconst_RMSE_curr = [resnormAll_lsqlin_unconst_RMSE_curr; output_inverse{curr_valToUse}(minRmseIndArray(curr_valToUse)).resnorm_lsqlin_unconst];
        end
        
        c_temp1 = mean(c_recovered_curr, 1);  c_temp1 = c_temp1/sum(c_temp1);
        J_temp1 = mean(J_recovered_curr, 1);  J_temp1 = J_temp1/sum(J_temp1);
        J_temp1_abs = mean(J_recovered_abs_curr, 1);
        lambda_temp1 = mean(lambda_recovered_curr, 1); 
        
        if sum(J_temp1) > 1.001
            akfdjalskf = 1;
        end
        
%         c_temp2 = mean(c_negConst_recovered_curr, 1);
        res_temp1 = mean(resnormAll_lsqlin_const_RMSE_curr, 1);
        rmse_temp1 = mean(rmseAll_lsqlin_const_RMSE_curr, 1);
%         res_temp2 = mean(resnormAll_lsqlin_unconst_RMSE_curr, 1);
        

%         c_temp2 = c_temp2/sum(c_temp2);
%         res_temp1 = res_temp1/sum(res_temp1);
%         res_temp2 = res_temp2/sum(res_temp2);
        
        avgWeightArray(ind_tCount, :) = c_temp1;
        avgRatioArray(ind_tCount, :) = J_temp1;
        avgAbsArray(ind_tCount, :) = J_temp1_abs;
        avgLambdaArray(ind_tCount, :) = lambda_temp1;
%         avgWeightNegConstArray(ind_tCount, :) = c_temp2;
        resnromAll_lsqlin_const_minRMSE_array(ind_tCount, :) = res_temp1;
        rmseAll_const_minRMSE_array(ind_tCount, :) = rmse_temp1; 
%         resnromAll_lsqlin_unconst_minRMSE_array(ind_tCount, :) = res_temp2;
        
        % remove the entries that are above the threshold
        keepInd = find(resnormAll_lsqlin_const_RMSE_curr < constThres);
        
        if isempty(keepInd)
            switch noWinMode
                case 'blank'
                    c_recover_belowThres = zeros(1, ccostCount);
                    J_recover_belowThres = zeros(1, ccostCount);
                    lambda_recover_belowThres = zeros(1, cconstCount);
                    resnorm_belowThres = 0;
                    rmse_belowThres = 0;
                    
%                 case 'previous'
%                     if ind_tCount == 1
%                         c_recover_belowThres = zeros(1, ccostCount);
%                         resnorm_belowThres = 0;
%                     else
%                         c_recover_belowThres = avgWeightArray_belowThres(ind_tCount - 1, :);
%                         resnorm_belowThres = resnromAll_lsqlin_const_minRMSE_array_belowThres(ind_tCount - 1, :);
%                     end
            end
        else
            c_temp = mean(c_recovered_curr(keepInd, :), 1);
            J_temp = mean(J_recovered_curr(keepInd, :), 1);
            lambda_temp = mean(lambda_recovered_curr(keepInd, :), 1);
            res_temp = mean(resnormAll_lsqlin_const_RMSE_curr(keepInd, :), 1);
            rmse_temp = mean(rmseAll_lsqlin_const_RMSE_curr(keepInd, :), 1);
            
            c_recover_belowThres = c_temp / (sum(c_temp));
            J_recover_belowThres = J_temp / (sum(J_temp));
            lambda_recover_belowThres = lambda_temp;
            resnorm_belowThres = res_temp;
            rmse_belowThres = rmse_temp;
        end
        
        avgWeightArray_belowThres(ind_tCount, :) = c_recover_belowThres;
        avgRatioArray_belowThres(ind_tCount, :) = J_recover_belowThres; 
        avgLambdaArray_belowThres(ind_tCount, :) = lambda_recover_belowThres;
        resnromAll_lsqlin_const_minRMSE_array_belowThres(ind_tCount, :) = resnorm_belowThres;
        resnromAll_lsqlin_const_minRMSE_array_belowThresPass(ind_tCount) = length(keepInd);
        rmse_const_minRMSE_array_belowThres(ind_tCount, :) = rmse_belowThres;
        
        % also pull out resnorm if it's only within seg bound
        if withinSegmentCheck(ind_tCount) 
             resnromAll_lsqlin_const_minRMSE_array_belowThres_inSegBound(ind_tCount, :) = resnorm_belowThres;
             rmse_const_minRMSE_array_belowThres_inSegBound(ind_tCount, :) = rmse_belowThres;
        else
             resnromAll_lsqlin_const_minRMSE_array_belowThres_inSegBound(ind_tCount, :) = 0;
             rmse_const_minRMSE_array_belowThres_inSegBound(ind_tCount, :) = 0;
        end
%         [motionEnd, motionWhole, ratioAll_end, ratioAll_all] = ...
%             checkStandingVSMovement(feature_full, t_windowCount, segmentInfo, ...
%             avgWeightArray_belowThres(ind_windowCount, :), ratioAll_end, ratioAll_all);
    else
        avgWeightArray(ind_tCount, :) = zeros(1, ccostCount);
        avgWeightNegConstArray(ind_tCount, :) = zeros(1, ccostCount);
        
        resnromAll_lsqlin_const_minRMSE_array(ind_tCount, :) = 0;
%         resnromAll_lsqlin_unconst_minRMSE_array(ind_tCount, :) = 0;
        
        avgWeightArray_belowThres(ind_tCount, :) = zeros(1, ccostCount);
        avgRatioArray_belowThres(ind_tCount, :) = zeros(1, ccostCount);
        avgLambdaArray_belowThres(ind_tCount, :) = zeros(1, cconstCount);
        resnromAll_lsqlin_const_minRMSE_array_belowThres(ind_tCount, :) = 0;
        resnromAll_lsqlin_const_minRMSE_array_belowThres_inSegBound(ind_tCount, :) = 0;
    end
    
     % now blend together the q over all the pass windows
     currT = feature_full.t(ind_tCount);     
     findTimeStep = find(t_belowThreshold_combined == currT);
     if ~isempty(findTimeStep) 
         % pull in existing values
         q_blended(:, ind_tCount) = mean(q_belowThreshold_combined(:, findTimeStep), 2);
         dq_blended(:, ind_tCount) = mean(dq_belowThreshold_combined(:, findTimeStep), 2);
         ddq_blended(:, ind_tCount) = mean(ddq_belowThreshold_combined(:, findTimeStep), 2);
         rmsEntryToUse(ind_tCount) = 1;
     elseif isempty(findTimeStep) && isempty(q_blended)
         % no value defined for that entry, and nothing to set as previous
         q_blended(:, ind_tCount) = zeros(size(q_belowThreshold_combined, 1), 1);
         dq_blended(:, ind_tCount) = zeros(size(dq_belowThreshold_combined, 1), 1);
         ddq_blended(:, ind_tCount) = zeros(size(ddq_belowThreshold_combined, 1), 1);
         rmsEntryToUse(ind_tCount) = 0;
     else
         switch noWinMode
             case 'blank'
                 q_blended(:, ind_tCount) = zeros(size(q_belowThreshold_combined, 1), 1);
                 dq_blended(:, ind_tCount) = zeros(size(dq_belowThreshold_combined, 1), 1);
                 ddq_blended(:, ind_tCount) = zeros(size(ddq_belowThreshold_combined, 1), 1);
                 rmsEntryToUse(ind_tCount) = 0;
                 
             case 'previous'
                 q_blended(:, ind_tCount) = q_blended(:, ind_tCount-1);
                 dq_blended(:, ind_tCount) = dq_blended(:, ind_tCount-1);
                 ddq_blended(:, ind_tCount) = ddq_blended(:, ind_tCount-1);
                 rmsEntryToUse(ind_tCount) = 1;
         end
     end
end

% indToUseStart = indToUse_window(1, 1); 
% indToUseEnd = indToUse_window(end, 2);
if sum(rmsEntryToUse) > 0
    blended_rmse_belowThreshold = rmse_fct(feature_full.q(:, find(rmsEntryToUse)), q_blended(:, find(rmsEntryToUse)), []);
    blended_rmse_belowThreshold_dq = rmse_fct(feature_full.dq(:, find(rmsEntryToUse)), dq_blended(:, find(rmsEntryToUse)), []);
    blended_rmse_belowThreshold_ddq = rmse_fct(feature_full.ddq(:, find(rmsEntryToUse)), ddq_blended(:, find(rmsEntryToUse)), []);
else
    blended_rmse_belowThreshold = -1;
    blended_rmse_belowThreshold_dq = -1;
    blended_rmse_belowThreshold_ddq = -1;
end

ccost_full_plot = ccost_full ./ repmat(sum(ccost_full, 1), size(ccost_full, 1), 1);

% resnorm rejection criteria
meanconstResnrom = mean(resnormAll_lsqlin_const_RMSE_plot);
% meanunconstResnorm = mean(resnormAll_lsqlin_unconst_RMSE_plot);

meanconstResnrom_averageWindow  = mean(resnromAll_lsqlin_const_minRMSE_array);
% meanunconstResnorm_averageWindow  = mean(resnormAll_lsqlin_unconst_RMSE_plot);

pts_aboveConstThres = resnormAll_lsqlin_const_RMSE_plot >= constThres;
pts_underConstThres = resnormAll_lsqlin_const_RMSE_plot <  constThres;

% pts_aboveDiffThres = resnormAll_lsqlin_const_RMSE_plot - resnormAll_lsqlin_unconst_RMSE_plot >= diffThreshold;
% pts_underDiffThres = resnormAll_lsqlin_const_RMSE_plot - resnormAll_lsqlin_unconst_RMSE_plot <  diffThreshold;

% pts_allow = or(pts_underConstThres, pts_underDiffThres);
% pts_deny = setxor(1:length(pts_allow), find(pts_allow));
% pts_deny =  pts_aboveConstThres && pts_aboveDiffThres;

pts_aboveConstThres_averageWindow = resnromAll_lsqlin_const_minRMSE_array >= constThres;
pts_underConstThres_averageWindow = resnromAll_lsqlin_const_minRMSE_array <  constThres;

% minWeightsUnconst = min(cAllUnconst_plot, [], 2);
% lowestValInUnconst = minWeightsUnconst > negConstThreshold;

% now we can do our metric counting
resnormPass1 = sum(pts_underConstThres);
resnormPass2 = length(resnormAll_lsqlin_const_RMSE_plot);
resnormPass3 = 0; % sum(lowestValInUnconst);
resnormPass4 = sum(pts_underConstThres_averageWindow);

weightsActive = sum(avgWeightArray_belowThres);
weightsMean = mean(avgWeightArray_belowThres);

% set them up in structs so we can tally them 
rmse_report.resnormPass1 = resnormPass1;
rmse_report.resnormPass2 = resnormPass2;
rmse_report.resnormPass3 = resnormPass3;
rmse_report.resnormPass4 = resnormPass4;

rmse_report.weightsActive = weightsActive;
rmse_report.weightsTotal = length(feature_full.t);
rmse_report.weightsMean = weightsMean;

rmse_report.resnorm_mean = mean(resnromAll_lsqlin_const_minRMSE_array_belowThres);
rmse_report.resnorm_std = std(resnromAll_lsqlin_const_minRMSE_array_belowThres);

resnromAll_lsqlin_const_minRMSE_array_belowThres_inSegBound_out = resnromAll_lsqlin_const_minRMSE_array_belowThres_inSegBound(resnromAll_lsqlin_const_minRMSE_array_belowThres_inSegBound > 0);
rmse_report.resnorm_mean_inseg = mean(resnromAll_lsqlin_const_minRMSE_array_belowThres_inSegBound_out);
rmse_report.resnorm_std_inseg = std(resnromAll_lsqlin_const_minRMSE_array_belowThres_inSegBound_out);

rmse_const_minRMSE_array_belowThres_inSegBound_out = rmse_const_minRMSE_array_belowThres_inSegBound(rmse_const_minRMSE_array_belowThres_inSegBound > 0);
rmse_report.rmse_mean_inseg = mean(rmse_const_minRMSE_array_belowThres_inSegBound_out);
rmse_report.rmse_std_inseg = std(rmse_const_minRMSE_array_belowThres_inSegBound_out);

rmse_report.blended_rmse_belowThreshold = blended_rmse_belowThreshold;
rmse_report.windowed_rmse_belowThreshold_mean = mean(windowed_rmse_belowThreshold);
rmse_report.windowed_rmse_belowThreshold_std = std(windowed_rmse_belowThreshold);

rmse_report.blended_rmse_belowThreshold_dq = blended_rmse_belowThreshold_dq;
rmse_report.windowed_rmse_belowThreshold_dq_mean = mean(windowed_rmse_belowThreshold_dq);
rmse_report.windowed_rmse_belowThreshold_dq_std = std(windowed_rmse_belowThreshold_dq);

rmse_report.blended_rmse_belowThreshold_ddq = blended_rmse_belowThreshold_ddq;
rmse_report.windowed_rmse_belowThreshold_ddq_mean = mean(windowed_rmse_belowThreshold_ddq);
rmse_report.windowed_rmse_belowThreshold_ddq_std = std(windowed_rmse_belowThreshold_ddq);