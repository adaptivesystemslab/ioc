function bool = checkSensorPlacement(mdl, sensorMarker, ...
    markerOffset, lengthUseInd, sensorMarkerStr, attachmentFrameStr, normThreshold)           

    if ~exist('normThreshold', 'var')
        normThreshold = 0.01;
    end

    bool = 1;
    vec = normVector(sensorMarker);
    vecStd = std(vec(lengthUseInd));
    
    if sum(vec(lengthUseInd)) == 0 || ... % marker missing in the entire lengthUseInd
            sum(sum(vec)) == 0 || ... % marker missing in the whole trajectory
            ~isempty(find(isnan(markerOffset))) || ... % NaN in the transformation matrix
            vecStd > normThreshold % it really shouldn't vary higher than 1 cm
        bool = 0;
        if ~isempty(sensorMarkerStr) && ~isempty(attachmentFrameStr)
            fprintf('WARNING: checkSensorPlacement was not able to attach sensor %s to frame %s due to sensor location malformation\n', ...
                sensorMarkerStr, attachmentFrameStr);
        end
    elseif ~isempty(mdl)
        try % run the same checks as 'obj.applyMarkerSensor'
            f = mdl.getFrameByName(attachmentFrameStr); % will crash if fails
        catch err
            fprintf('WARNING: checkSensorPlacement did not find frame: %s\n', attachmentFrameStr);
            bool = 0;
        end
    end