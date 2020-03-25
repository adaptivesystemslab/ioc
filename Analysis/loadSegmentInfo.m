function [segData, segOnlyDataTable, restOnlyDataTable] = loadSegmentInfo(filepathSegments, trialInfo)
% filepathSegments = 'D:\aslab\data\Fullbody_IIT_2017\ManualSeg.xlsx';
    manSegTab = readtable(filepathSegments,'ReadRowNames',true);
    
    startColStr = ['S' trialInfo.runName(8:9) '_START'];
    endColStr = ['S' trialInfo.runName(8:9) '_END'];
    
    for i = 1:14
        segData(i).state = manSegTab.State(i);
        segData(i).direction = manSegTab.Direction(i);
        segData(i).count = manSegTab.Count(i);
        
        if ~isempty(manSegTab.(startColStr)(i))
            segData(i).timeStart = manSegTab.(startColStr)(i);
            segData(i).timeEnd = manSegTab.(endColStr)(i);
        else
            segData(i).timeStart = 0;
            segData(i).timeEnd = 0;
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