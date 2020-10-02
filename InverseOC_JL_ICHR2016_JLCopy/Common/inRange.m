function bool = inRange(valInsp, valTarget, errorRange)
    % bool = inRange(valInsp, valTarget, errorRange)
    %
    % determines if 'valInsp' is within 'errorRange' or 'valTarget'
    
    if ~exist('errorRange', 'var')
        % seems small enough
        errorRange = 1e-3;
    end

    valTemp = abs(valInsp - valTarget);
    
    if valTemp <= errorRange
        bool = 1;
    else
        bool = 0;
    end
end