function numGradSum = calcNumGradSum(depFctsUncheck, indepVarsUncheck, indepValsUncheck, h)
    if ~iscell(depFctsUncheck)
        depFcts{1} = depFctsUncheck;
    else
        depFcts = depFctUncheck;
    end

    if ~iscell(indepVarsUncheck)
        indepVars{1} = indepVarsUncheck;
    else
        indepVars = indepVarsUncheck;
    end
    
    if ~iscell(indepValsUncheck)
        indepVals = num2cell(indepValsUncheck, size(indepValsUncheck, 2));
    else
        indepVals = indepValsUncheck;
    end

    for ind_preallocate = 1:length(indepVars{1})
        numGradSum(ind_preallocate) = sym(0);
    end
    
    for ind_baseFct = 1:length(depFcts) % iterate through all the dependent variables
        for ind_indepVar = 1:length(indepVars) % iterate through all the independent variables
            currNumJac = calcNumGrad(depFcts{ind_baseFct}, indepVars{ind_indepVar}, indepVals{ind_indepVar}, h);
            numGradSum = numGradSum + currNumJac; % sum all the different gradient functions together
        end % ind_indepVar
    end % ind_baseFct
end

function numGrad = calcNumGrad(depFcts, indepVars, indepVals, h)
    % use the finite diffrence formula to calculate the numerical gradient
    for ind_preallocate = 1:length(indepVars)
        numGrad(ind_preallocate) = sym(0);
    end
    
    for ind_indepVars = 1:length(indepVars)
        currIndepVar = indepVars(ind_indepVars); % select one of the variables to differentiate
        currIndepVal = indepVals(ind_indepVars);
        
        [missIndepVar, IA, ~] = setxor(indepVars, currIndepVar); % and build an array that excludes the selected variable
        missIndepVal = indepVals(IA);
        
        indep_ph = currIndepVal + h; % calculate the value of the current varible +/- h
        indep_mh = currIndepVal - h;
        dep_ph = subs(depFcts, currIndepVar, indep_ph); % as well as the corresponding value in the function
        dep_mh = subs(depFcts, currIndepVar, indep_mh);
        
        numGradSym = (dep_ph - dep_mh) / (indep_ph - indep_mh);
        numGrad(ind_indepVars) = subs(numGradSym, missIndepVar, missIndepVal);
    end
end