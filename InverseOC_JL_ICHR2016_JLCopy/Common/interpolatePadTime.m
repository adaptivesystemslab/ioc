function outData = interpolatePadTime(targetTime, sourceTime, sourceData)
    % take in target time, and interpolate source time where interpolation
    % is possible. for extrapolation components, the sourceData is padded 
    % instead
    
    interpolationmode = 'interpolateandpad'; % spline interpolateandpad
    
    padVal = 0; 
    initTime = min([targetTime(1) sourceTime(1)]);
    dof = size(sourceData, 2);
    
    targetTime = targetTime - initTime;
    sourceTime = sourceTime - initTime;
    
    if 0
        figure;
        plot([targetTime(1) targetTime(end)], [0 0], 'r');
        hold on
        plot([sourceTime(1) sourceTime(end)], [1 1], 'b');
        ylim([-0.5 1.5]);
    end
    
    if isequal(targetTime, sourceTime)
        outData = sourceData;
        return
    end
    
        x = sourceTime';
    y = sourceData';
    xx = targetTime; 
    sourceSpline = spline(x, y, xx)';
    
        % check the dims...
    if dof ~= size(sourceSpline, 2)
        sourceSpline = sourceSpline';
    end
    
    switch interpolationmode
        case 'spline'
            
        case 'interpolateandpad'
            [closeValStart, closeIndStart] = findClosestValue(sourceTime(1), targetTime);
            [closeValEnd, closeIndEnd] = findClosestValue(sourceTime(end), targetTime);
            
            if closeIndStart > 1
                sourceSpline(1:closeIndStart, 1:dof) = padVal;
            end
            
            if closeIndEnd < length(targetTime)
                xtopad = length(targetTime) - closeIndEnd;
                sourceSpline(closeIndEnd:length(targetTime), 1:dof) = padVal;
            end
    end
    
    outData = sourceSpline;
    
%     find(sourceTime > targetTime(1));
%     
%     % from target(1) to start of source(1) - padding
%     if sourceTime(1) < targetTime(1)
%         % the source array is larger than the target - crop
%         
%     elseif sourceTime(1) == targetTime(1)
%         % source array is equal size - do nothing for now
%         
%     else
%         % the source array is smaller than the target - padding
%     end
%     
%     % from source(1) to source(end) - interpolate
%     
%     % from source(end) to target(end) - padding
%     