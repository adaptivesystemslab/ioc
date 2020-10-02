function    pointsToRemove = findStartEndPoints(segLabel, objLabel)
    % produces an array the contains the values for the start and end
    searchInd = find(segLabel == objLabel);
    diffSearchInd = [1; diff(searchInd)];
    gapSearchArrayMid = find(diffSearchInd > 1);

    % figure out the points to remove
    if length(gapSearchArrayMid) > 0
        pointsToRemoveStart = 1:gapSearchArrayMid(1);
    
%     pointsToRemoveEnd = searchInd(gapSearchArrayMid(end)):searchInd(end);
        pointsToRemoveEnd = [];
    
    	pointsToRemove = [pointsToRemoveStart pointsToRemoveEnd];
    else
        pointsToRemove = [];
    end