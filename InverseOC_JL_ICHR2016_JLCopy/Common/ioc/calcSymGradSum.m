function symGradSum = calcSymGradSum(depFctsUncheck, indepVarsUncheck, indepVals)
    % perform symbolic gradient calculations
    % depFcts [cell array]: the dependent functions
    % indepVars [cell array]: the independent variables
    
    if ~iscell(depFctsUncheck)
        depFcts{1} = depFctsUncheck;
    else
        depFcts = depFctsUncheck;
    end
    
    if ~iscell(indepVarsUncheck)
        indepVars{1} = indepVarsUncheck;
    else
        indepVars = indepVarsUncheck;
    end
    
    for ind_preallocate = 1:length(indepVars{1})
        symGradSum(ind_preallocate) = sym(0);
    end

    for ind_baseFct = 1:length(depFcts) % iterate through all the dependent variables
        for ind_indepVar = 1:length(indepVars) % iterate through all the independent variables
            currSymJac = calcSymGrad(depFcts{ind_baseFct}, indepVars{ind_indepVar});
            symGradSum = symGradSum + currSymJac; % sum all the different gradient functions together
        end % ind_indepVar
    end % ind_baseFct
    
    if ~isempty(indepVals)
        for ind_indepVar = 1:length(indepVars) % iterate through all the independent variables
            symGradSum = subs(symGradSum, indepVars{ind_indepVar}, indepVals{ind_indepVar}); % sub in the indep values
        end % ind_indepVar
    end
end

function symGrad = calcSymGrad(depFcts, indepVars)
    % use the MATLAB function 'jacobian' to calculate the symbolic gradient
    symGrad = jacobian(depFcts, indepVars);
end