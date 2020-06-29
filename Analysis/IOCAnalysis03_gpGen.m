%restructure the mat file for GP analysis
%% load the mat in question
loadFilePath = 'D:\results\fatigue_ioc03_weightsPattern\20200413_FatigueFull_3CF_3\mat\mat_Subject02_3DOF_3CF.mat';
load(loadFilePath);

%% rearange into x = rep, y = features
% pull out the q information
segmentInfoTable = struct2table(segmentInfo);
segmentInfoTable.Properties.VariableNames{'state'} = 'segType';
[segOnlyDataTable, restOnlyDataTable, segMask, restMask] = sepSegRest(segmentInfoTable);
fatigueScore = segOnlyDataTable.fatigueScore;

y = [];
y_label = {};
y_inds = 5:22;
for i = 1:3
    segStats_SingleWindow = stats_q(i).segStats_SingleWindow;
    [segOnlyDataTable, restOnlyDataTable, segMask, restMask] = sepSegRest(segStats_SingleWindow);
    
    x = 1:size(segOnlyDataTable, 1);
    y_fieldNames = segOnlyDataTable.Properties.VariableNames;
    y_seg = table2array(segOnlyDataTable(:, y_inds));
    y = [y y_seg];
    for j = y_inds
        y_new_label = [stats_q(i).name '_' y_fieldNames{j}];
        y_label = [y_label y_new_label];
    end
end

 
y = [y fatigueScore];
y_label = [y_label 'fatigueScore'];