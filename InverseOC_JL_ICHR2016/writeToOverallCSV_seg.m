function writeToOverallCSV(matFilesToLoad, exportPath, currDataset, currInstName, elapsedTime, runSettings, rmse_report, cost_function_names, param, comments, thresholdMultiplier, constThres, weightCounter, metricsToMeasure) 
    if ~exist(exportPath, 'file')
        % if doesn't exist, write header
        header = ['runFile,dataset,instance,type,comments,thresholdMultiplier,constThres,win_length,win_shift,spline_length,'...
            'win_length_sim,doc_sim_win_len_sim,spline_length_sim,', ...
            'half_phase_construction,doc_sim_sp_dq_endcond,doc_sim_cf_dq_constraint,doc_sim_win_length,ioc_knots,doc_knots,correlation_threshold,' ...
            'ioc_add_knot_as_cf_constraints,ioc_sp_dq_endcond,ioc_cf_q_constraint,ioc_cf_dq_constraint,ioc_knot_locations,segment_only_windows,segStart,segEnd,accMax,dofMax,thresMax,dirMax,tpMax,tnMax,fnMax,fpMax'];
        
        for i = 1:length(cost_function_names)
            header = [header ',acc_' cost_function_names{i}];
        end
        
                for i = 1:length(cost_function_names)
            header = [header ',maxsel_' cost_function_names{i}];
        end
        
        header = [header ',meanRMSE,rangeRMSE,stdRMSE,runtime,lenTime,resnormpass1sum,totalnumberofresnormpasspossible,resormpass3sum,resormpass4sum,', ...
            'rmse_blended,rmse_windowed_mean,rmse_windowed_std,resnorm_mean,resnorm_std,resnorm_mean_seg,resnorm_std_seg,rmse_windowed_mean_seg,rmse_windows_std_seg'];
        
         if ~isempty(metricsToMeasure)
        for i = 1:length(metricsToMeasure)
            header = [header ',' metricsToMeasure{i}.name];
        end
         end
        
    else
        header = '';
    end
    
    fId = fopenCheck(exportPath, header);
    
    fprintf(fId, ['%s,%s,%s,%s,%s,%f,%f,%u,%u,%u,' ...
        '%u,%f,%u,' ...
        '%s,%s,%s,%u,%u,%u,%f,' ...
        '%s,%s,%s,%s,%s,%s,%f,%f'], ...
        matFilesToLoad{1}, currDataset, currInstName, currInstName(end-12:end-1), comments, thresholdMultiplier, constThres, param.win_length, param.win_shift, param.spline_length, ...
        runSettings.variableFactors.win_length_sim, runSettings.variableFactors.doc_sim_win_length_sim, runSettings.variableFactors.spline_length_sim, ...
        param.half_phase_construction, param.doc_sim_sp_dq_endcond, param.doc_sim_cf_dq_constraint, param.doc_sim_win_length, param.ioc.knot_count, param.doc_pivot.knot_count, param.corrThreshold, ...
        param.ioc_add_knot_as_cf_constraints, param.ioc_sp_dq_endcond, param.ioc_cf_q_constraint, param.ioc_cf_dq_constraint, param.ioc_knot_locations, param.segment_only_windows, rmse_report.clusterRange(1), rmse_report.clusterRange(2));
    
    fprintf(fId, ',%s', cost_function_names{rmse_report.dofMax});
    fprintf(fId, ',%f', rmse_report.accMax(rmse_report.dofMax).acc);
    fprintf(fId, ',%f', rmse_report.thresMax(rmse_report.dofMax));
    fprintf(fId, ',%s', rmse_report.dirMax{rmse_report.dofMax});
    
    fprintf(fId, ',%f', rmse_report.accMax(rmse_report.dofMax).tp); % x,tpMax,tnMax,fnMax,fpMax'];
    fprintf(fId, ',%f', rmse_report.accMax(rmse_report.dofMax).tn);
    fprintf(fId, ',%f', rmse_report.accMax(rmse_report.dofMax).fn);
    fprintf(fId, ',%f', rmse_report.accMax(rmse_report.dofMax).fp);
    
    for i = 1:length(rmse_report.accMax)
        outputStr = [',%f'];
        fprintf(fId, outputStr, rmse_report.accMax(i).acc);
    end
    
    for i = 1:length(rmse_report.accMax)
        outputStr = [',%f'];
        fprintf(fId, outputStr, weightCounter(i));
    end
    
    fprintf(fId, [',%f,%f,%f,%f,%f,%f,%f,%f,%f,' ...
        '%f,%f,%f,%f,%f,' ...
        '%f,%f,%f,%f'], ...
        rmse_report.meanRMSE, rmse_report.rangeRMSE, rmse_report.stdRMSE, 0, rmse_report.t(end), rmse_report.resnormPass1, rmse_report.resnormPass2, rmse_report.resnormPass3, rmse_report.resnormPass4, ...
        rmse_report.blended_rmse_belowThreshold, rmse_report.windowed_rmse_belowThreshold_mean, rmse_report.windowed_rmse_belowThreshold_std, rmse_report.resnorm_mean, rmse_report.resnorm_std, ...
        rmse_report.resnorm_mean_inseg, rmse_report.resnorm_std_inseg, rmse_report.rmse_mean_inseg, rmse_report.rmse_std_inseg);
    
    if ~isempty(metricsToMeasure)
    for i = 1:length(metricsToMeasure)
        outputStr = [',%f'];
        fprintf(fId, outputStr, metricsToMeasure{i}.value);
    end
    end
    
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