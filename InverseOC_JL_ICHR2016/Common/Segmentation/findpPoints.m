function pPoints = findpPoints(tLabel, targetVal, pointLimit)
    pPointsInit = find(tLabel == targetVal);

    if pointLimit > length(pPointsInit)
        pointLimit = length(pPointsInit);
    end

    pDs = sort(randperm(length(pPointsInit), pointLimit)); % p downsample
    pPoints = pPointsInit(pDs);
end