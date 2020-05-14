function normVec = normVector(origMatrix, param)
    % normVec = normVector(origMatrix); 
    % calculate the distance vector of a given matrix, into a vector
    
    if ~exist('param', 'var')
        param.allowNeg = 0;
    end
    
    len = size(origMatrix, 1);
    
    if param.allowNeg
        % subtract the min so it allows negative values from the norm, thus
        % removes 'rebound' effects when it goes below zero
        normInit = norm(origMatrix(1, :)); % will reset the norm so it's back to this as starting value
        origMatrix = origMatrix - min(min(origMatrix));
    end
    
    normVec = zeros(len, 1);
    for x = 1:len
%         sq = origMatrix(x, :) .^ 2;
%         normVec(x) = sqrt(sum(sq));

        normVec(x) = norm(origMatrix(x, :));
    end
    
    if param.allowNeg
        % restore the offset applied previously by shifting the initial
        % time series point to where it would've been if it wasn't for the
        % offset
        normVec = normVec - normVec(1) + normInit; 
    end