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
C_report{2+len_name*0, 1} = '-c_orig-';
C_report{3+len_name*1, 1} = '-c-';
C_report{4+len_name*2, 1} = '-J-';
C_report{5+len_name*3, 1} = '-cJ-';

for ii_report = 1:len_name
    C_report{1, 3+ii_report} = ['C_' cost_function_names{ii_report}];
    C_report{2+ii_report, 1} = cost_function_names{ii_report};
    C_report{3+len_name+ii_report, 1} = cost_function_names{ii_report};
    C_report{4+len_name*2+ii_report, 1} = cost_function_names{ii_report};
    C_report{5+len_name*3+ii_report, 1} = cost_function_names{ii_report};
end

C_report(2+len_name*0+(1:len_name), 2) = mat2cell(ccost_win', ones(length(ccost_win), 1), 1);
C_report(3+len_name*1+(1:len_name), 2) = mat2cell(ccost_win'/sum(ccost_win), ones(length(ccost_win), 1), 1);
C_report(4+len_name*2+(1:len_name), 2) = mat2cell(doc_optsim_win.J_array', ones(length(doc_optsim_win.J_array), 1), 1);
C_report(5+len_name*3+(1:len_name), 2) = mat2cell(doc_optsim_win.J_contrib', ones(length(doc_optsim_win.J_contrib), 1), 1);

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
    
    currCOriginal = output_inverse{ind_windowCount}(ii_touse).c_original;
    currCRecovered = output_inverse{ind_windowCount}(ii_touse).c_recovered;
    currJVal = output_inverse{ind_windowCount}(1).J_out_array;
    currContribVal = J_inverse{ii_touse};
    
    C_report(2+len_name*0+(1:len_name), ii_report+2) = mat2cell(currCOriginal, 1, ones(1, len_name));
    C_report(3+len_name*1+(1:len_name), ii_report+2) = mat2cell(currCRecovered, 1, ones(1, len_name));
    C_report(4+len_name*2+(1:len_name), ii_report+2) = mat2cell(currJVal, 1, ones(1, len_name));
    C_report(5+len_name*3+(1:len_name), ii_report+2) = mat2cell(currContribVal, 1, ones(1, len_name));
end

ind_sectionPre = 6+len_name*4-1;

C_report{ind_sectionPre+1, 1} = '-metrics-';
C_report{ind_sectionPre+2, 1} = 'rmse(q)';
C_report{ind_sectionPre+3, 1} = 'norm(res)';
C_report{ind_sectionPre+4, 1} = 'norm(res) orig';
C_report{ind_sectionPre+5, 1} = 'resnormScaleFactor';

for ii_report = 1:length(ccost_win)+1
    if ii_report == 1
        ii_touse = minRmseInd; % pull out the entry with the smallest RMSE
    else
        ii_touse = ii_report - 1; % everything else is normal as expected
    end
    
    currRMSE = output_inverse{ind_windowCount}(ii_touse).rmse; 
    curr_resnorm_lsqlin = output_inverse{ind_windowCount}(ii_touse).resnorm_lsqlin; 
    curr_resnorm_lsqlin_orig = output_inverse{ind_windowCount}(ii_touse).resnorm_orig;
    curr_condNum = output_inverse{ind_windowCount}(ii_touse).condA;
    curr_resnormScale = output_inverse{ind_windowCount}(ii_touse).c_scalingFactor;
    
    C_report(ind_sectionPre+2, ii_report+2) = mat2cell(currRMSE, ones(length(currRMSE), 1), 1);
    C_report(ind_sectionPre+3, ii_report+2) = mat2cell(curr_resnorm_lsqlin, ones(length(curr_resnorm_lsqlin), 1), 1);
    C_report(ind_sectionPre+4, ii_report+2) = mat2cell(curr_resnorm_lsqlin_orig, ones(length(curr_condNum), 1), 1);
    C_report(ind_sectionPre+5, ii_report+2) = mat2cell(curr_resnormScale, ones(length(curr_condNum), 1), 1);
end

% C_report{ind_sectionPre+5, 1} = 'lambda(A)';
% C_report(9+len_name*3:8+len_name*3+len_lambda, end-len_name:end) = ...
%     mat2cell([lambda_array(:, minRmseInd) lambda_array], ones(1, len_lambda), ones(1, len_name+1));