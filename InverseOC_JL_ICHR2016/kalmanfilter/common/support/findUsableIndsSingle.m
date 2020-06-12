function usableInds = findUsableIndsSingle(correctArray, markerStr, modelInstance, featureSet)
    markerInd = find(ismember(featureSet.measurement_labels, markerStr));
    usableInds = find(correctArray(:, markerInd));
end