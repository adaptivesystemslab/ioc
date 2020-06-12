function [corrCosts, rref_out, eig_vals, eig_cutoff] = ...
    checkCorrelationSVD(J, param, J_out_contrib_array_full, functionMode)
    % want to ensure that the selected basis functions are uniquely
    % determined
    if ~exist('functionMode', 'var') || isempty(functionMode)
        functionMode = 'prevWin';
    end
    
    functionModeIfEmpty = 'largest'; 
    functionIfLargest = 'cost'; % grad, cost
    len_cost_functions = size(J, 2);
    
    svd_ratio_threshold = param.corrThreshold;
    corrCosts_prev = param.cost_functions_ioc;
    
    [rref_out, jb, eig_cutoff, eig_vals] = calcSVDThreshold(J, 1:len_cost_functions, svd_ratio_threshold, len_cost_functions, param);
    defaultArray = 1:len_cost_functions;
%     max_rref_out = max(rref_out);
%     abssum_rref_out = sum(abs(rref_out));
    count_rref_out = sum(abs(rref_out) > 0);
%     [findrow_rref_out, findcol_rref_out] = find(rref_out == 1);

    if 0
        % produce the array that has 1 at valid cost functions and 0 at dep
        % ones, which we'll estimate by taking the mean of the matrix
        corrCosts = defaultArray(max_rref_out >= 1);
    else
       % produce the array that aims to minimize the overall J magnitude
%        linearDepCol = find(max_rref_out < 1); 
        linearDepCol = find(count_rref_out > 1);
%         linearDepCol = [find(count_rref_out == 0) find(count_rref_out > 1)];
        
        %        assessJ = sum(abs(J));
        switch functionIfLargest
            case 'grad'
                assessJ = normVector(J')';
                 
            case 'cost'
                assessJ = J_out_contrib_array_full;
        end
        
        % if there is a 0 in count_rref_out, it is dependent and needs to
        % be dropped (TODO)
        toDropInd = find(count_rref_out == 0);
        defaultArray(toDropInd) = 0;
        assessJ(toDropInd) = 0; %
        
        
       switch functionMode
           case 'largest'
               corrCosts = modeHighest(linearDepCol, rref_out, assessJ, defaultArray);
               
           case 'prevWin'              
               if isempty(corrCosts_prev) || length(corrCosts_prev) == len_cost_functions
                   % if previous unique functions were not defined, or were 
                   % all independent, then just go with an existing option
                   corrCosts = modeHighest(linearDepCol, rref_out, assessJ, defaultArray);
                   
               else
                   % otherwise, take the previous window
                   for ind = 1:length(linearDepCol)
                       linearDepCol_curr = rref_out(:, linearDepCol(ind));
                       overlapEntries = [find(abs(linearDepCol_curr) > 0); linearDepCol(ind)]';
                       
                       % so from within overlapEntries, one of the items needs to be
                       % dropped. identify all the dep entries that were
                       % not indep from previous windows
                       toDropInd = setdiff(overlapEntries, corrCosts_prev);
                       
                       if isempty(toDropInd)
                           % no entries is getting dropped, then go with
                           % the maxJ idea
                           [~, toDropInd] = max(assessJ(overlapEntries));
                       end
                       
                       defaultArray(toDropInd) = 0;
                       assessJ(toDropInd) = 0; % once it's been dropped, can't pick it again
                   end
                   
                    corrCosts = defaultArray(defaultArray > 0);
               end
       end
    end
end

function [rref_out, jb, eig_cutoff, eig_vals] = calcSVDThreshold(J, defaultArray, svd_ratio_threshold, len_cost_functions, param)

    J_test = J(:, defaultArray);
    maxJ = max(abs(J_test));
    [sortVal, sortInd] = sort(maxJ); %  B = A(IX)
    [sortValRev, ~] = sort(maxJ, 'descend');
    
    if svd_ratio_threshold > 0
        [eig_cutoff, eig_vals, eig_perc] = screePlotMethod(J_test, svd_ratio_threshold);
    else
        [~, eig_vals, eig_perc] = screePlotMethod(J_test, 1);
        eig_cutoff = len_cost_functions + svd_ratio_threshold*1;
    end
    
    if svd_ratio_threshold >= 1 || eig_cutoff >= len_cost_functions 
        % don't cut out any values
        [rref_out, jb] = rref(J_test, eps);
    else
        % cut out values by threshold. pick the next smallest cutoff and
        % scale it up a little since putting the cutoff right at the
        % selected value may cause numerical roundoff problems
        
        cutoffTol = mean([sortValRev(eig_cutoff) sortValRev(eig_cutoff+1)]); % use the maxabs because that's what rref is using
        [rref_out1, jb1] = rref(J_test(:, sortInd), cutoffTol);
        [rref_out2, jb2] = rref(J_test(:, fliplr(sortInd)), cutoffTol);
        [rref_out3, jb3] = rref(J_test(:, :), cutoffTol);
        
        name1 = param.cost_function_names(sortInd);
        name2 = param.cost_function_names(fliplr(sortInd));
        name3 = param.cost_function_names;
        
        % restore the original order
        rref_out = zeros(size(rref_out1));
        for ind = 1:length(sortInd)
            rref_out(:, sortInd(ind)) = rref_out1(:, ind);
            name(sortInd(ind)) = name1(ind);
        end
        
        jb = sortInd(jb1);
    end
end

function corrCosts = modeHighest(linearDepCol, rref_out, assessJ, defaultArray)
    for ind = 1:length(linearDepCol)
        linearDepCol_curr = rref_out(:, linearDepCol(ind));
        overlapEntries = [find(abs(linearDepCol_curr) > 0); linearDepCol(ind)];

        % so from within overlapEntries, one of the items needs to be
        % dropped. pick the largest J
        [~, toDropInd] = max(assessJ(overlapEntries));
        defaultArray(toDropInd) = 0;
        assessJ(toDropInd) = 0; % once it's been dropped, can't pick it again
    end

    corrCosts = defaultArray(defaultArray > 0);
end