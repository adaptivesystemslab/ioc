function [segOnlyDataTable, restOnlyDataTable, segMask, restMask] = sepSegRest(tableOrig)
    segMask = find(ismember(tableOrig.segType, 'Seg'));
    segOnlyDataTable = tableOrig(segMask, :);
    
    restMask = find(ismember(tableOrig.segType, 'Rest'));
    restOnlyDataTable = tableOrig(restMask, :);
end