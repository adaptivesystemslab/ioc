function [closeValArray, closeIndArray, targetArray, diffVal] = findClosestValue(valueToFindInput, targetArray, favour)
    % [closeVal, closeInd] = findClosestValue(valueToFind, targetArray, favour)
    %
    % finds the closest value. 'favour' is an optional variable
    % - 'above' will return the closest value larger than valueToFind
    % - 'below' will return the closest value smaller than valueToFind
    
    % automatically sorts the array
    targetArray = sort(targetArray);
    closeValArray = zeros(size(valueToFindInput));
    closeIndArray = zeros(size(valueToFindInput));
    
    for i = 1:length(valueToFindInput)
        valueToFind = valueToFindInput(i);
        closestOver = find(targetArray - valueToFind > 0, 1, 'first');
        closestExact = find(targetArray == valueToFind, 1); % very possible!
        closestUnder = find(targetArray - valueToFind < 0, 1, 'last');
        %     closestUnder = closestOver - 1;
        
        if ~isempty(closestExact)
            % exact match in array
            closeInd = closestExact;
            
        elseif ~isempty(closestOver) && isempty(closestUnder)
            % only top side matched
            closeInd = closestOver;
            
        elseif isempty(closestOver) && ~isempty(closestUnder)
            % only bottom side matched
            closeInd = closestUnder;
            
        elseif isempty(closestOver) && isempty(closestUnder)
            % no match
            closeInd = length(targetArray);
            
        else
            % there is one match in both sides
            if exist('favour', 'var')
                % does user have a direction preference?
                switch favour
                    case 'above'
                        closeInd = closestOver;
                        
                    case 'below'
                        closeInd = closestUnder;
                        
                    otherwise
                        % defaulting to below, if jibberish
                        closeInd = closestUnder;
                end
                
            else
                % take the actual closest one
                overMag = abs(targetArray(closestOver) - valueToFind);
                underMag = abs(targetArray(closestUnder) - valueToFind);
                
                if overMag > underMag
                    closeInd = closestUnder;
                else
                    closeInd = closestOver;
                end
            end
        end
        
        closeVal = targetArray(closeInd);
        
        closeValArray(i) = closeVal;
        closeIndArray(i) = closeInd;        
    end
    
    diffVal = valueToFindInput - closeValArray;
end