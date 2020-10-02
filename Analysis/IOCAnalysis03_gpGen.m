%restructure the mat file for GP analysis
%% load the mat in question
loadFilePath = 'C:\Users\jf2lin\Downloads\TransferXL-00YNMdn2vnv70\png\mat\mat_Subject14_3DOF_3CF.mat';
load(loadFilePath);

setPaths();

%% rearange into x = rep, y = features
% pull out the q information
segmentInfoTable = struct2table(segmentInfo);
segmentInfoTable.Properties.VariableNames{'state'} = 'segType';
[segOnlyDataTable, restOnlyDataTable, segMask, restMask] = sepSegRest(segmentInfoTable);
fatigueScore = segOnlyDataTable.fatigueScore;

y = [];
y_label = {};
y_inds = 5:22;
% for i = 1:3
%     currStats = stats_q(i);
%     segStats_SingleWindow = currStats.segStats_SingleWindow;
%     [segOnlyDataTable, restOnlyDataTable, segMask, restMask] = sepSegRest(segStats_SingleWindow);
%     
%     x = 1:size(segOnlyDataTable, 1);
%     y_fieldNames = segOnlyDataTable.Properties.VariableNames;
%     y_seg = table2array(segOnlyDataTable(:, y_inds));
%     y = [y y_seg];
%     for j = y_inds
%         y_new_label = [currStats.name '_' y_fieldNames{j}];
%         y_label = [y_label y_new_label];
%     end
% end
% 
% for i = 1:3
%     currStats = stats_dq(i);
%     segStats_SingleWindow = currStats.segStats_SingleWindow;
%     [segOnlyDataTable, restOnlyDataTable, segMask, restMask] = sepSegRest(segStats_SingleWindow);
%     
%     x = 1:size(segOnlyDataTable, 1);
%     y_fieldNames = segOnlyDataTable.Properties.VariableNames;
%     y_seg = table2array(segOnlyDataTable(:, y_inds));
%     y = [y y_seg];
%     for j = y_inds
%         y_new_label = [currStats.name '_' y_fieldNames{j}];
%         y_label = [y_label y_new_label];
%     end
% end
% 
% for i = 1:3
%     currStats = stats_tau(i);
%     segStats_SingleWindow = currStats.segStats_SingleWindow;
%     [segOnlyDataTable, restOnlyDataTable, segMask, restMask] = sepSegRest(segStats_SingleWindow);
%     
%     x = 1:size(segOnlyDataTable, 1);
%     y_fieldNames = segOnlyDataTable.Properties.VariableNames;
%     y_seg = table2array(segOnlyDataTable(:, y_inds));
%     y = [y y_seg];
%     for j = y_inds
%         y_new_label = [currStats.name '_' y_fieldNames{j}];
%         y_label = [y_label y_new_label];
%     end
% end
% 
% for i = 1:3
%     currStats = stats_dtau(i);
%     segStats_SingleWindow = currStats.segStats_SingleWindow;
%     [segOnlyDataTable, restOnlyDataTable, segMask, restMask] = sepSegRest(segStats_SingleWindow);
%     
%     x = 1:size(segOnlyDataTable, 1);
%     y_fieldNames = segOnlyDataTable.Properties.VariableNames;
%     y_seg = table2array(segOnlyDataTable(:, y_inds));
%     y = [y y_seg];
%     for j = y_inds
%         y_new_label = [currStats.name '_' y_fieldNames{j}];
%         y_label = [y_label y_new_label];
%     end
% end

for i = 1:3
    currStats = stats_weights(i);
    segStats_SingleWindow = currStats.segStats_SingleWindow;
    [segOnlyDataTable, restOnlyDataTable, segMask, restMask] = sepSegRest(segStats_SingleWindow);
    
    x = 1:size(segOnlyDataTable, 1);
    y_fieldNames = segOnlyDataTable.Properties.VariableNames;
    y_seg = table2array(segOnlyDataTable(:, y_inds));
    y = [y y_seg];
    for j = y_inds
        y_new_label = [currStats.name '_' y_fieldNames{j}];
        y_label = [y_label y_new_label];
    end
end

for i = 1:3
    currStats = stats_dweights(i);
    segStats_SingleWindow = currStats.segStats_SingleWindow;
    [segOnlyDataTable, restOnlyDataTable, segMask, restMask] = sepSegRest(segStats_SingleWindow);
    
    x = 1:size(segOnlyDataTable, 1);
    y_fieldNames = segOnlyDataTable.Properties.VariableNames;
    y_seg = table2array(segOnlyDataTable(:, y_inds));
    y = [y y_seg];
    for j = y_inds
        y_new_label = [currStats.name '_' y_fieldNames{j}];
        y_label = [y_label y_new_label];
    end
end

y = [y fatigueScore];
y_label = [y_label 'fatigueScore'];

saveVar.x = x;
saveVar.y = y;
saveVar.y_label = y_label;

% save('saveVar.mat', 'saveVar');