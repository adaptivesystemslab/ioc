function [currExerType, currExerInst] = exerciseStringDelimit(currExer)
    % given the exercise string, figure out which part is the name of the
    % exercise and which part is the instance value
    
    currExerSplit = strsplit(currExer, '_');
    intCheck = str2num(currExerSplit{end});
    if ~isempty(intCheck)
        % remove the final underscore and note the instance
        currExerType = currExer(1:end-1-length(currExerSplit{end}));
        currExerInst = intCheck;
    else
        % otherwise, there's no underscore at the end. iterate
        % through the name array to extract the name and the
        % instance
        for ind_exerciseNameLength = 1:length(currExer)
            % want to remove the trailing number on the exercise name
            intCheck = str2num(currExer(ind_exerciseNameLength:end));
            if ~isempty(intCheck)
                currExerType = currExer(1:ind_exerciseNameLength-1);
                currExerInst = str2num(currExer(ind_exerciseNameLength:end));
                break
            end
        end
    end
