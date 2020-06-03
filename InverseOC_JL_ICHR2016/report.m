len_name = length(cost_function_names);
len_ran = length(param.cost_functions_ioc);
len_svd = length([output_inverse{ind_windowCount}(1).svd_s]);
len_rmse = length(currRMSE_set(1).arrayName);
len_lambda = size(lambda_array, 1);
len_Jbreakdown = length(direct_check_flags{1}.J_breakdown.report);
len_residual = length(output_inverse{ind_windowCount}(1).residual_lsqlin);

C_report = {};
% C_report{1, 1} = ['t=' num2str(t_opt(1)) '-'  num2str(t_opt(end))];

C_report{1, 2} = 'C_gen';
C_report{1, 3} = ['C_' cost_function_names{minRmseInd} '_' num2str(indToUse_window(ind_windowCount, 1)) '_'  num2str(indToUse_window(ind_windowCount, 2))];
C_report{2, 1} = '-cost-';
C_report{3+len_name, 1} = '-contrib-';
for ii_report = 1:len_name
    C_report{1, 3+ii_report} = ['C_' cost_function_names{ii_report}];
    C_report{2+ii_report, 1} = cost_function_names{ii_report};
    C_report{3+len_name+ii_report, 1} = cost_function_names{ii_report};
end

C_report(2+(1:len_name), 2) = mat2cell(ccost_win', ones(length(ccost_win), 1), 1);
C_report(3+len_name+(1:len_name), 2) = mat2cell(doc_optsim_win.J_contrib', ones(length(doc_optsim_win.J_contrib), 1), 1);

for ii_report = 1:length(ccost_win)+1
    if ii_report == 1
        ii_touse = minRmseInd; % pull out the entry with the smallest RMSE
    else
        ii_touse = ii_report - 1; % everything else is normal as expected
    end
    
    % if this entry was not ran due to dependences, fill this section with
    % blanks
    [closeVal, closeInd] = findClosestValue(ii_touse, param.cost_functions_ioc);
    
    if closeVal ~= ii_touse
        % we don't have a match
        continue
    end
    
    ccost_offset = ccost_win(ii_touse);
    if ccost_offset == 0
        ccost_offset = 1;
    end

    currPivotVal = ccost_inverse_pivot{ii_touse}*ccost_offset;
    currContribVal = J_inverse{ii_touse};
    C_report(2+(1:len_name), ii_report+2) = mat2cell(currPivotVal, 1, ones(1, len_name));
    C_report(3+len_name+(1:len_name), ii_report+2) = mat2cell(currContribVal, 1, ones(1, len_name));
end

C_report{4+len_name*2, 1} = '-metrics-';

% C_report{5+len_name*2, 1} = 'rmse(q_opt)';
% C_report(5+len_name*2, end-len_name:end) = mat2cell([rmse(minRmseInd) rmse], 1, ones(1, len_name+1));

C_report(5+len_name*2:5+len_name*2+len_rmse-1, 1) = currRMSE_set(1).arrayName;
C_report(5+len_name*2:5+len_name*2+len_rmse-1, end-len_name:end) = mat2cell([rmse_array(:, minRmseInd) rmse_array], ones(1, len_rmse), ones(1, len_name+1));
C_report{5+len_name*2+len_rmse, 1} = 'count(rmse)';
C_report(5+len_name*2+len_rmse, end-len_name:end) = mat2cell([minIndCount(minRmseInd) minIndCount], 1, ones(1, len_name+1));

C_report{6+len_name*2+len_rmse, 1} = 'norm(res)';
resnorm_refit = [output_inverse{ind_windowCount}(:).resnorm];
C_report(6+len_name*2+len_rmse, end-len_name:end) = ...
    mat2cell([output_inverse{ind_windowCount}(minRmseInd).resnorm resnorm_refit], 1, ones(1, len_name+1));

C_report{7+len_name*2+len_rmse, 1} = 'DOC_exitflag';
C_report(7+len_name*2+len_rmse, 2) = mat2cell(doc_optsim_win.exitflag', ones(length(doc_optsim_win.exitflag), 1), 1);

C_report{8+len_name*2+len_rmse, 1} = 'IOC_exitflag';

C_report{9+len_name*2+len_rmse, 1} = 'J';
C_report(9+len_name*2+len_rmse, 2) = mat2cell(doc_optsim_win.J', ones(length(doc_optsim_win.J), 1), 1);

C_report{10+len_name*2+len_rmse, 1} = 'c';
C_report(10+len_name*2+len_rmse, 2) = mat2cell(doc_optsim_win.c', ones(length(doc_optsim_win.c), 1), 1);

C_report{11+len_name*2+len_rmse, 1} = 'ceq_x';
C_report(11+len_name*2+len_rmse, 2) = mat2cell(doc_optsim_win.ceq_x', ones(length(doc_optsim_win.ceq_x), 1), 1);

C_report{12+len_name*2+len_rmse, 1} = 'ceq_dx';
C_report(12+len_name*2+len_rmse, 2) = mat2cell(doc_optsim_win.ceq_dx', ones(length(doc_optsim_win.ceq_dx), 1), 1);

C_report{13+len_name*2+len_rmse, 1} = 'ceq_ddx';
C_report(13+len_name*2+len_rmse, 2) = mat2cell(doc_optsim_win.ceq_ddx', ones(length(doc_optsim_win.ceq_ddx), 1), 1);

% C_report{21+len_svd+len_name*2+len_rmse+len_lambda-2, 1} = 'Jbreakdown+Ind';
% C_report((21+len_svd+len_name*2+len_rmse+len_lambda-2):(21+len_svd+len_name*2+len_rmse+len_lambda-2+len_Jbreakdown-1), 2) = mat2cell(doc_optsim_win.J_breakdown.report, ones(length(doc_optsim_win.J_breakdown.report), 1), 1);

C_report{(21+len_svd+len_name*2+len_rmse+len_lambda-2+len_Jbreakdown), 1} = 'resnorm_lsqlin';
C_report{(21+len_svd+len_name*2+len_rmse+len_lambda-2+len_Jbreakdown+1), 1} = 'residual_lsqlin';


for ii_report = 1:length(ccost_win)+1
    if ii_report == 1
        ii_touse = minRmseInd; % pull out the entry with the smallest RMSE
    else
        ii_touse = ii_report - 1; % everything else is normal as expected
    end
    
    currDOCExitFlag = direct_check_flags{ii_touse}.exitflag;
    currIOCExitFlag = output_inverse{ind_windowCount}(ii_touse).exitflag; 
    curr_resnorm_lsqlin = output_inverse{ind_windowCount}(ii_touse).resnorm_lsqlin; 
    curr_residual_lsqlin = output_inverse{ind_windowCount}(ii_touse).residual_lsqlin; 
    currJ = direct_check_flags{ii_touse}.J;
    currC = direct_check_flags{ii_touse}.c;
    currCeqX = direct_check_flags{ii_touse}.ceq_x;
    currCeqDx = direct_check_flags{ii_touse}.ceq_dx;
    currCeqDdx = direct_check_flags{ii_touse}.ceq_ddx;
    currJ_breakdown = direct_check_flags{ii_touse}.J_breakdown.report;
    
    C_report(7+len_name*2+len_rmse, ii_report+2) = mat2cell(currDOCExitFlag, ones(length(currDOCExitFlag), 1), 1);
    C_report(8+len_name*2+len_rmse, ii_report+2) = mat2cell(currIOCExitFlag, ones(length(currIOCExitFlag), 1), 1);
    C_report(9+len_name*2+len_rmse, ii_report+2) = mat2cell(currJ, ones(length(currJ), 1), 1);
    C_report(10+len_name*2+len_rmse, ii_report+2) = mat2cell(currC, ones(length(currC), 1), 1);
    C_report(11+len_name*2+len_rmse, ii_report+2) = mat2cell(currCeqX, ones(length(currCeqX), 1), 1);
    C_report(12+len_name*2+len_rmse, ii_report+2) = mat2cell(currCeqDx, ones(length(currCeqDx), 1), 1);
    C_report(13+len_name*2+len_rmse, ii_report+2) = mat2cell(currCeqDdx, ones(length(currCeqDdx), 1), 1);

    currGradDotSign = output_inverse{ind_windowCount}(ii_touse).gradDotSign;
    currSvd = output_inverse{ind_windowCount}(ii_touse).svd_s;
    C_report{14+len_name*2+len_rmse, ii_report+2} = num2str(currGradDotSign);
%     C_report((19+len_name*2+len_rmse):(19+len_svd+len_name*2-1+len_rmse), ii_report+2) = mat2cell(currSvd, ones(1, len_svd), ones(1, 1));

%     C_report((21+len_svd+len_name*2+len_rmse+len_lambda-2):(21+len_svd+len_name*2+len_rmse+len_lambda-2+len_Jbreakdown-1), ii_report+2) = mat2cell(currJ_breakdown, ones(length(currJ_breakdown), 1), 1);
% 
    C_report((21+len_svd+len_name*2+len_rmse+len_lambda-2+len_Jbreakdown), ii_report+2) = mat2cell(curr_resnorm_lsqlin, ones(length(curr_resnorm_lsqlin), 1), 1);
%     C_report((21+len_svd+len_name*2+len_rmse+len_lambda-2+len_Jbreakdown+1):(21+len_svd+len_name*2+len_rmse+len_lambda-2+len_Jbreakdown+len_residual), ii_report+2) ...
%         = mat2cell(curr_residual_lsqlin, ones(length(curr_residual_lsqlin), 1), 1);
end

C_report{14+len_name*2+len_rmse, 1} = 'sign(pivotdot)';

C_report{15+len_name*2+len_rmse, 1} = '-cond-';

C_report{16+len_name*2+len_rmse, 1} = 'cond(A)';
cond_refit = [output_inverse{ind_windowCount}(:).condA];
C_report(16+len_name*2+len_rmse, end-len_name:end) = ...
        mat2cell([output_inverse{ind_windowCount}(minRmseInd).condA cond_refit], 1, ones(1, len_name+1));
    
C_report{17+len_name*2+len_rmse, 1} = 'corr(A)';

% pull out the top correlation values that will fit into the space

C_report{19+len_name*2+len_rmse, 1} = 'svd(A)';

C_report{20+len_svd+len_name*2-1+len_rmse, 1} = 'lambda';
C_report(20+len_svd+len_name*2-1+len_rmse:20+len_svd+len_name*2-1+len_rmse+len_lambda-1, end-len_name:end) = ...
    mat2cell([lambda_array(:, minRmseInd) lambda_array], ones(1, len_lambda), ones(1, len_name+1));

% C_report{22, 1} = 'J_res';
% C_report{22, 2} = res_direct;
% C_report(22, end-length(ccost_inverse_pivot)+1:end) = mat2cell(res_inverse, 1, ones(1, length(ccost_inverse_pivot)));

% % settings report
% C_metrics = {};
% C_metrics{1, 2} = 'C_gen';
% C_metrics{1, 3} = 'C_pivot_ddq';
% C_metrics{1, 4} = 'C_pivot_tau';
% C_metrics{2, 1} = '-DC exit-';
% C_metrics{3, 1} = '-IOC exit-';
% C_metrics{4, 1} = '-DC exit-';
% C_metrics{2, 2} = direct_opt_flags(1).Exitflag;
% for ii_blah = 1:length(ccost_inverse_pivot)
%     C_metrics{3, ii_blah+2} = inverse_flags{ii_blah}.Exitflag;
%     C_metrics{4, ii_blah+2} = direct_check_flags{ii_blah}.Exitflag;
% end

% C_report(:, 1) = [];
C_report_array{ind_windowCount} = C_report;
% C_metrics_array{ind_windowCount} = C_metrics;

