function [clusterStart, clusterEnd, clusterLength] = findDiscontinuity(array, mergeFactor)
% given an array, look for breaks in the continuity

% % short way
%         diffSearchInd = [1; diff(array)];
%     gapSearchArrayMid = find(diffSearchInd > 1);
    
    % long way
    prevVal = array(1);
    clusterCount = 1;
    clusterStart = 1;
    
    for ind = 2:length(array)
        currVal = array(ind);
        if sum(prevVal + (1:mergeFactor) == currVal)
            % do nothing
        else
            % declare a new cluster
            clusterEnd(clusterCount) = ind-1;
            
            clusterCount = clusterCount + 1;
            clusterStart(clusterCount) = ind;
        end
        prevVal = currVal;
    end
    
    clusterEnd(clusterCount) = length(array);
    clusterLength = clusterEnd - clusterStart;
    
