function [segData, segOnlyDataTable, restOnlyDataTable] = loadSegmentInfo(filepathSegments, trialInfo)
% filepathSegments = 'D:\aslab\data\Fullbody_IIT_2017\ManualSeg.xlsx';
    manSegTab = readtable(filepathSegments,'ReadRowNames',true);
    
    startColStr = ['S' trialInfo.runName(8:9) '_START'];
    endColStr = ['S' trialInfo.runName(8:9) '_END'];
    scoreColStr = ['S' trialInfo.runName(8:9) '_SCORE'];
    
    % look for the starting and ending rest
    restCount = 0;
    
    tableInd = 0;
    for i = 1:size(manSegTab, 1)
            tableInd = tableInd + 1;
            segData(tableInd).state = manSegTab.State(i);
%             segData(tableInd).direction = manSegTab.Direction(i);
            segData(tableInd).count = manSegTab.Count(i);
            segData(tableInd).timeStart = manSegTab.(startColStr)(i);
            segData(tableInd).timeEnd = manSegTab.(endColStr)(i);
            segData(tableInd).fatigueScore = manSegTab.(scoreColStr)(i);
%             
            if strcmpi(manSegTab.State(i), 'rest') && segData(tableInd).fatigueScore == 8
                   break;
            end
            
%             if isnan(manSegTab.(scoreColStr)(i))
%                 segData(tableInd).fatigueScore = manSegTab.(scoreColStr)(i-1);
%             else
%                 segData(tableInd).fatigueScore = manSegTab.(scoreColStr)(i);
%             end
    end
    
    for i = 1:length(segData)
        if strcmpi(segData(i).state, 'rest') ...
                && isnan(segData(i).timeStart) && isnan(segData(i).timeEnd)
            segData(i).timeStart = segData(i-1).timeEnd;
            segData(i).timeEnd = segData(i+1).timeStart;
            segData(i).fatigueScore = segData(i-1).fatigueScore;
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