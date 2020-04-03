function [segData, segOnlyDataTable, restOnlyDataTable] = loadSegmentInfo(filepathSegments, trialInfo)
% filepathSegments = 'D:\aslab\data\Fullbody_IIT_2017\ManualSeg.xlsx';
    manSegTab = readtable(filepathSegments,'ReadRowNames',true);
    
    startColStr = ['S' trialInfo.runName(8:9) '_START'];
    endColStr = ['S' trialInfo.runName(8:9) '_END'];
    scoreColStr = ['S' trialInfo.runName(8:9) '_SCORE'];
    
    tableInd = 0;
    for i = 1:size(manSegTab, 1)
        if ~isnan(manSegTab.(startColStr)(i))
            tableInd = tableInd + 1;
            segData(tableInd).state = manSegTab.State(i);
%             segData(tableInd).direction = manSegTab.Direction(i);
            segData(tableInd).count = manSegTab.Count(i);
            segData(tableInd).timeStart = manSegTab.(startColStr)(i);
            segData(tableInd).timeEnd = manSegTab.(endColStr)(i);
            
            if isnan(manSegTab.(scoreColStr)(i))
                segData(tableInd).fatigueScore = manSegTab.(scoreColStr)(i-1);
            else
                segData(tableInd).fatigueScore = manSegTab.(scoreColStr)(i);
            end
            
        end
    end
    
%     figure;
%     plot(trajT, q);
%     title(savePath);
 
    segDataTable = struct2table(segData);
    mask = ismember(segDataTable{:,'state'}, 'Seg');
    segOnlyDataTable = segData(mask);
%     plotBoxes(segOnlyDataTable, [0 1 0], 0.1);
    
    mask = ismember(segDataTable{:,'state'}, 'Rest');
    restOnlyDataTable = segData(mask);
%     plotBoxes(segOnlyDataTable, [1 0 0], 0.1);
    

end