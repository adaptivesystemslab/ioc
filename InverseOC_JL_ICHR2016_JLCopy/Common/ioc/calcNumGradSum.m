function [numGradSum, numGradSym] = calcNumGradSum_sym(depFctsUncheck, indepVarsUncheck, indepValsUncheck, h)
    syms h_sym real
    
    if ~iscell(depFctsUncheck)
        depFcts{1} = depFctsUncheck;
    else
        depFcts = depFctsUncheck; 
    end

    if ~iscell(indepVarsUncheck)
        indepVars{1} = indepVarsUncheck; sf
    else
        indepVars = indepVarsUncheck;
    end
    
    if isempty(indepValsUncheck)
        % empty. we will return a symbolic expression
        indepVals{1} = [];        
    elseif ~iscell(indepValsUncheck)
        indepVals{1} = indepValsUncheck;
    else
        indepVals = indepValsUncheck;
    end

    numGradSym = sym(zeros(size(indepVars{1})));
    numGradSum = sym(zeros(size(indepVals{1})));
    
    for ind_baseFct = 1:length(depFcts) % iterate through all the dependent variables
        for ind_indepVar = 1:length(indepVars) % iterate through all the independent variables
            currNumJac_Sym = ...
                calcNumGrad(depFcts{ind_baseFct}, indepVars{ind_indepVar}, h_sym);
            
            numGradSym = numGradSym + currNumJac_Sym;
            
            % if passed in a numerical array, we'll evaluate it
            if ~isempty(numGradSum)
                localIndepVars = [indepVars{ind_indepVar} h_sym];
                
                for ind_time = 1:size(indepVals{ind_indepVar}, 1)
                    localIndepVals = [indepVals{ind_indepVar}(ind_time, :) h]; % TODO iterate through time
                    currNumJac = subs(currNumJac_Sym, localIndepVars, localIndepVals);
                    
                    numGradSum(ind_time, :) = numGradSum(ind_time, :) + currNumJac; % sum all the different gradient functions together
                end
            end
        end % ind_indepVar
    end % ind_baseFct
end

function numGradSym = calcNumGrad(depFcts, indepVars, h)
%     % use the finite diffrence formula to calculate the numerical gradient
%     for ind_preallocate = 1:length(indepVars)
%         for ind_len = 1:size(indepVals, 1)
%             numGrad(ind_preallocate) = sym(0);
%         end
%     end
    
    for ind_indepVars = 1:length(indepVars)
        currIndepVar = indepVars(ind_indepVars); % select one of the variables to differentiate
        
        indep_ph = currIndepVar + h; % calculate the value of the current varible +/- h
        indep_mh = currIndepVar - h;
        dep_ph = subs(depFcts, currIndepVar, indep_ph); % as well as the corresponding value in the function
        dep_mh = subs(depFcts, currIndepVar, indep_mh);
        
        numGradSym(ind_indepVars) = (dep_ph - dep_mh) ./ (indep_ph - indep_mh);
    end
end