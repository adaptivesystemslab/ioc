function [closeVal, closeInd] = findClosestValueInplace(valueToFind, targetArray, favour)
    % [closeVal, closeInd] = findClosestValue(valueToFind, targetArray, favour)
    %
    % finds the closest value. 'favour' is an optional variable
    % - 'above' will return the closest value larger than valueToFind
    % - 'below' will return the closest value smaller than valueToFind
    
    % automatically sorts the array
%     targetArray = sort(targetArray);
    
    closestOver = find(targetArray - valueToFind > 0);
    closestExact = find(targetArray == valueToFind, 1); % very possible!
    closestUnder = find(targetArray - valueToFind < 0);
%     closestUnder = closestOver - 1;

    if ~isempty(closestExact)
        % exact match in array
        closeInd = closestExact;
        
    elseif ~isempty(closestOver) && isempty(closestUnder)
        % only top side matched
        if length(closestOver) > 1
            % multiple matches, since it's over, take the smallest
            targetArrayNew = targetArray(closestOver);
            [y, i] = min(targetArrayNew);
            closeInd = find(targetArray == y);
        else
            % only one matched
            closeInd = closestOver;
        end
        
    elseif isempty(closestOver) && ~isempty(closestUnder)
        % only bottom side matched
        if length(closestUnder) > 1
            % multiple matches, since it's under, take the largest
            targetArrayNew = targetArray(closestUnder);
            [y, i] = max(targetArrayNew);
            closeInd = find(targetArray == y);
        else
            closeInd = closestUnder;
        end
        
    elseif isempty(closestOver) && isempty(closestUnder)
        % no match
        closeInd = length(targetArray);
        
    else
%         if exist('favour', 'var')
%             % does user have a direction preference?
%             switch favour
%                 case 'above'
%                     closeInd = closestOver;
%                     
%                 case 'below'
%                     closeInd = closestUnder;
%                     
%                 otherwise
%                     % defaulting to below
%                     closeInd = closestUnder;
%             end
%             
%         else
            overMag = abs(targetArray(closestOver) - valueToFind);
            underMag = abs(targetArray(closestUnder) - valueToFind);

            if overMag > underMag
                closeInd = closestUnder;
            else
                closeInd = closestOver;
            end
%         end
    end
    
    closeVal = targetArray(closeInd);
end