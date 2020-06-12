function writeToInstanceCSV(exportPath, cost_function_names, startTime, endTime, ...
    set_ioc_check_struct, ccost_win, doc_optsim_win, ccost_inverse_pivot, ...
    J_inverse, currRMSE_set, output_inverse, ind_windowCount, minRmseInd, direct_check_flags, motionWhole, param) 
    if ~exist(exportPath, 'file') 
        % if doesn't exist, write header  
        header = ['runTime,t_start,t_end,qCk_Rmse,qDoc_Rmse']; % indToUse_window(ind_windowCount, 1), indToUse_window(ind_windowCount, 2)
        
        for i = 1:length(cost_function_names)
            header = [header ',C_gen_cost_' cost_function_names{i}];
        end
        
        for i = 1:length(cost_function_names)
            header = [header ',C_gen_const_' cost_function_names{i}];
        end
        
        header = [header ',optpivot'];
        
        for i = 1:length(cost_function_names)
            header = [header ',C_optpivot_cost_' cost_function_names{i}];
        end
        
        for i = 1:length(cost_function_names)
            header = [header ',C_optpivot_const_' cost_function_names{i}];
        end
        
%         for i = 1:length(cost_function_names)
%             header = [header ',C_optpivot_costunconst_' cost_function_names{i}];
%         end
        
        header = [header ',rmse, resnorm, IOC_exitflag, J_DOCRMSE, c, ceq_x, ceq_dx, ceq_ddx, condA, motionWholeProfile'];
    else
        header = '';
    end
    
    fId = fopenCheck(exportPath, header);
    ii_touse = minRmseInd;
      
%   header = ['runTime,t_start,t_end,qCk_Rmse,qDoc_Rmse,']; % indToUse_window(ind_windowCount, 1), indToUse_window(ind_windowCount, 2)
    outputStr = ['%s,%f,%f,%s,%s'];
    fprintf(fId, outputStr, datestr(now), startTime, endTime, ...
        arrayToStr(set_ioc_check_struct(ind_windowCount).intermed_ind_set + startTime-1), ...
        arrayToStr(set_ioc_check_struct(ind_windowCount).const_x + startTime-1));
    
    for i = 1:length(cost_function_names)
%         header = [header ',C_gen_cost_' cost_function_names{i}];
        outputStr = [',%f'];
        fprintf(fId, outputStr, ccost_win(i));
    end
    
    for i = 1:length(cost_function_names)
%         header = [header ',C_gen_const_' cost_function_names{i}];
        outputStr = [',%f'];
        fprintf(fId, outputStr, doc_optsim_win.J_contrib(i));
    end
    
%     header = [header ',optpivot'];
    outputStr = [',%s'];
    fprintf(fId, outputStr, cost_function_names{minRmseInd});
    
    currPivotVal = ccost_inverse_pivot{minRmseInd};
    currContribVal = J_inverse{minRmseInd};
    
    for i = 1:length(cost_function_names)
%         header = [header ',C_optpivot_cost_' cost_function_names{i}];
        outputStr = [',%f'];
        fprintf(fId, outputStr, currPivotVal(i));
    end
    
    for i = 1:length(cost_function_names)
%         header = [header ',unconst' cost_function_names{i}];
        outputStr = [',%f'];
        fprintf(fId, outputStr, currContribVal(i));
    end
    
%     c_out_pivot_unconst = output_inverse{ind_windowCount}(ii_touse).c_out_pivot_unconst;
%     
%     for i = 1:length(cost_function_names)
%         outputStr = [',%f'];
%         fprintf(fId, outputStr, c_out_pivot_unconst(i));
%     end
        
%     header = [header ',rmse, resnorm, IOC_exitflag, J_DOCRMSE, c, ceq_x, ceq_dx, ceq_ddx, condA, \n'];
  
    currDOCExitFlag = direct_check_flags{ii_touse}.exitflag;

    currIOCExitFlag = output_inverse{ind_windowCount}(ii_touse).exitflag;
    curr_resnorm_lsqlin = output_inverse{ind_windowCount}(ii_touse).resnorm_lsqlin;
%     curr_resnorm_unconst_lsqlin = output_inverse{ind_windowCount}(ii_touse).resnorm_lsqlin_unconst;
%     curr_residual_lsqlin = output_inverse{ind_windowCount}(ii_touse).residual_lsqlin;
    currJ = direct_check_flags{ii_touse}.J;
    currC = direct_check_flags{ii_touse}.c;
    currCeqX = direct_check_flags{ii_touse}.ceq_x;
    currCeqDx = direct_check_flags{ii_touse}.ceq_dx;
    currCeqDdx = direct_check_flags{ii_touse}.ceq_ddx;
    condA = output_inverse{ind_windowCount}(ii_touse).condA;
%     ioc_corrJHqmax_overall = output_inverse{ind_windowCount}(ii_touse).corrJHqmax_overall;
%     ioc_corrJHqmin_overall = output_inverse{ind_windowCount}(ii_touse).corrJHqmin_overall; 
%     ioc_corrJHabsmax_overall = output_inverse{ind_windowCount}(ii_touse).corrJHabsmax_overall;
%     ioc_corrJHqrange_overall = output_inverse{ind_windowCount}(ii_touse).corrJHqrange_overall;
%     
%     ioc_corrJHqmax_indiv = output_inverse{ind_windowCount}(ii_touse).corrJHqmax_indiv;
%     ioc_corrJHqmin_indiv = output_inverse{ind_windowCount}(ii_touse).corrJHqmin_indiv;
%     ioc_corrJHabsmax_indiv = output_inverse{ind_windowCount}(ii_touse).corrJHabsmax_indiv;
%     ioc_corrJHqrange_indiv = output_inverse{ind_windowCount}(ii_touse).corrJHqrange_indiv;
            
    outputStr = [',%f,%f,%f,%f,%f,%f,%f,%f,%f,%s \n'];
    fprintf(fId, outputStr, currRMSE_set(minRmseInd).q, curr_resnorm_lsqlin, currIOCExitFlag, currJ, ...
        currC, currCeqX, currCeqDx, currCeqDdx, condA, motionWhole);
    
    fclose(fId);
end