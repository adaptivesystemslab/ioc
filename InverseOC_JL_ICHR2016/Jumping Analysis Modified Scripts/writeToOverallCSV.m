function writeToOverallCSV(exportPath, loadPath, currInstName, elapsedTime, runSettings, rmse_report, cost_function_names, param)
    cost_function_names_resnorm = [cost_function_names, 'resnorm'];
    cost_function_names_lambda = rmse_report.rsq_namel;

    if ~exist(exportPath, 'file')
        % if doesn't exist, write header
        header = ['runTime,loadPath,instance,type,win_length,win_shift,spline_length,'...
            'win_length_sim,doc_sim_win_len_sim,spline_length_sim,', ...
            'half_phase_construction,doc_sim_sp_dq_endcond,doc_sim_cf_dq_constraint,doc_sim_win_length,ioc_knots,doc_knots,correlation_threshold,' ...
            'ioc_add_knot_as_cf_constraints,ioc_sp_dq_endcond,ioc_cf_q_constraint,ioc_cf_dq_constraint,ioc_knot_locations,segment_only_windows,' ...
            'meanRMSE,rangeRMSE,stdRMSE,runtime,lenTime,resnormpass1sum,resnormpass2sum,resormpass3sum,resormpass4sum,', ...
            'rmse_blended,rmse_windowed_mean,rmse_windowed_std,resnorm_mean,resnorm_std,resnorm_mean_seg,resnorm_std_seg,rmse_windowed_mean_seg,rmse_windows_std_seg'];
        
        for i = 1:length(cost_function_names)
            header = [header ',weightssum_' cost_function_names{i}];
        end
        
        for i = 1:length(cost_function_names)
            header = [header ',weightsmean_' cost_function_names{i}];
        end
        
        for i = 1:length(cost_function_names_resnorm)
            header = [header ',r_mean_fatigue_' cost_function_names_resnorm{i}];
        end
        
        for i = 1:length(cost_function_names_resnorm)
            header = [header ',r_std_fatigue_' cost_function_names_resnorm{i}];
        end
        
        for i = 1:length(cost_function_names_resnorm)
            header = [header ',r_mean_set5_' cost_function_names_resnorm{i}];
        end
        
        for i = 1:length(cost_function_names_resnorm)
            header = [header ',r_std_set5_' cost_function_names_resnorm{i}];
        end

        
%         for i = 1:length(cost_function_names_lambda)
%             header = [header ',r_mean_fatiguel_' cost_function_names_lambda{i}];
%         end
%         
%         for i = 1:length(cost_function_names_lambda)
%             header = [header ',r_std_fatiguel_' cost_function_names_lambda{i}];
%         end
%         
%         for i = 1:length(cost_function_names_lambda)
%             header = [header ',r_mean_set5l_' cost_function_names_lambda{i}];
%         end
%         
%         for i = 1:length(cost_function_names_lambda)
%             header = [header ',r_std_set5l_' cost_function_names_lambda{i}];
%         end

        
        header = [header ',total'];
        
    else
        header = '';
    end
    
    fId = fopenCheck(exportPath, header);
    
    fprintf(fId, ['%s,%s,%s,%s,%u,%u,%u,' ...
        '%u,%f,%u,' ...
        '%s,%s,%s,%u,%u,%u,%f,' ...
        '%s,%s,%s,%s,%s,%s,' ...
        '%f,%f,%f,%f,%f,%f,%f,%f,%f,' ...
        '%f,%f,%f,%f,%f,' ...
        '%f,%f,%f,%f'], ...
        runSettings.nowTimeStr, loadPath, currInstName, currInstName, param.win_length, param.win_shift, param.spline_length, ...
        runSettings.variableFactors.win_length_sim, runSettings.variableFactors.doc_sim_win_length_sim, runSettings.variableFactors.spline_length_sim, ...
        param.half_phase_construction, param.doc_sim_sp_dq_endcond, param.doc_sim_cf_dq_constraint, param.doc_sim_win_length, param.ioc.knot_count, param.doc_pivot.knot_count, param.corrThreshold, ...
        param.ioc_add_knot_as_cf_constraints, param.ioc_sp_dq_endcond, param.ioc_cf_q_constraint, param.ioc_cf_dq_constraint, param.ioc_knot_locations, param.segment_only_windows, ...
        rmse_report.meanRMSE, rmse_report.rangeRMSE, rmse_report.stdRMSE, elapsedTime, rmse_report.t(end), rmse_report.resnormPass1, rmse_report.resnormPass2, rmse_report.resnormPass3, rmse_report.resnormPass4, ...
        rmse_report.blended_rmse_belowThreshold, rmse_report.windowed_rmse_belowThreshold_mean, rmse_report.windowed_rmse_belowThreshold_std, rmse_report.resnorm_mean, rmse_report.resnorm_std, ...
        rmse_report.resnorm_mean_inseg, rmse_report.resnorm_std_inseg, rmse_report.rmse_mean_inseg, rmse_report.rmse_std_inseg);
    
    for i = 1:length(cost_function_names)
        outputStr = [',%f'];
        fprintf(fId, outputStr, rmse_report.weightsActive(i));
    end
    
    for i = 1:length(cost_function_names)
        outputStr = [',%f'];
        fprintf(fId, outputStr, rmse_report.weightsMean(i));
    end
    
    for i = 1:length(cost_function_names_resnorm)
        outputStr = [',%f'];
        fprintf(fId, outputStr, rmse_report.rsq_mean_fatigue(i));
    end
    
    for i = 1:length(cost_function_names_resnorm)
        outputStr = [',%f'];
        fprintf(fId, outputStr, rmse_report.rsq_std_fatigue(i));
    end
    
    for i = 1:length(cost_function_names_resnorm)
        outputStr = [',%f'];
        fprintf(fId, outputStr, rmse_report.rsq_mean_set5(i));
    end
    
    for i = 1:length(cost_function_names_resnorm)
        outputStr = [',%f'];
        fprintf(fId, outputStr, rmse_report.rsq_std_set5(i));
    end
    
    
%     for i = 1:length(cost_function_names_lambda)
%         outputStr = [',%f'];
%         fprintf(fId, outputStr, rmse_report.rsq_mean_fatiguel(i));
%     end
%     
%     for i = 1:length(cost_function_names_lambda)
%         outputStr = [',%f'];
%         fprintf(fId, outputStr, rmse_report.rsq_std_fatiguel(i));
%     end
%     
%     for i = 1:length(cost_function_names_lambda)
%         outputStr = [',%f'];
%         fprintf(fId, outputStr, rmse_report.rsq_mean_set5l(i));
%     end
%     
%     for i = 1:length(cost_function_names_lambda)
%         outputStr = [',%f'];
%         fprintf(fId, outputStr, rmse_report.rsq_std_set5l(i));
%     end
    
    
    fprintf(fId, ',%f', rmse_report.weightsTotal);
    
    fprintf(fId, '\n');
    
    fclose(fId);
end

% set them up in structs so we can tally them 
% rmse_report.resnormPass1 = resnormPass1;
% rmse_report.resnormPass2 = resnormPass2;
% rmse_report.resnormPass3 = resnormPass3;
% 
% rmse_report.weightsActive = weightsActive;
% rmse_report.weightsMean = weightsMean;