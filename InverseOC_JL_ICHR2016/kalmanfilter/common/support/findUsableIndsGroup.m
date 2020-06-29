function usableInds = findUsableIndsGroup(correctArray, jointCentreName, modelInstance, featureSet)
    % find that joint centre's composition markers
    jointCentreInd = [];
    for i = 1:length(modelInstance.jointCentreArray)
        if strcmpi(modelInstance.jointCentreArray{i}{1}, jointCentreName)
            jointCentreInd = i;
            break;
        end
    end
    
    masterUseIndMatrix = [];
    numMarkers = length(modelInstance.jointCentreArray{jointCentreInd}) - 1;
    for j = 2:length(modelInstance.jointCentreArray{jointCentreInd})
        markerStr = modelInstance.jointCentreArray{jointCentreInd}{j};
        
        switch markerStr
            case {'HARRINGTON_R', 'HARRINGTON_L'}
                markerData = ones(length(featureSet.time), 1);
                masterUseIndMatrix = [masterUseIndMatrix markerData];
                
            otherwise
                markerInd = find(ismember(featureSet.measurement_labels, markerStr));
                masterUseIndMatrix = [masterUseIndMatrix correctArray(:, markerInd)];
        end
    end

    sumCol = sum(masterUseIndMatrix, 2);
    usableInds = find(sumCol == numMarkers);
end