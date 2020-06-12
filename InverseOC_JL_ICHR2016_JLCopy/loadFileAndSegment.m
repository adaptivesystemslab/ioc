function [totalIndWeighting, accLess, accMore, cost_function_names, constThres, maxInd, rmse_report] = ...
    loadFileAndSegment(matFilesToLoad, currInstName, indTotal, outputPathLocal, masterOverall, masterSaveOverallSeg, comments, thresholdMultiplier, masterFigureSave)

plotFig = 0;

    outputPath = outputPathLocal;
    checkMkdir(outputPath);
    rmse_fct = @(a,b,c) sqrt(sum(sum((a-b).^2))/(size(a, 1)*size(a, 2)));
    
    nowStr = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
    nowStr = '';
    
    dof = 9;
    
    commentSplit = strsplit(comments, '_');
    offset = str2num(commentSplit{1});
    mansegMode = commentSplit{2}; % fullseg halfseg

    plotCalc_assembleFile;
    plotCalc_metricCalc;
    
    rmse_report.t = feature_full.t;
    rmse_report.maxminRMSE = [];
    rmse_report.meanRMSE = [];
    rmse_report.rangeRMSE = [];
    rmse_report.stdRMSE = [];
    
    switch currDataset
        case 'squats_tuat_2011'
            switch mansegMode
                case 'halfseg'
                    % if it's fatigue dataset, we want to create intermed segments
                    shiftZVCSettings.offsetVal = -0.0; % initial offset applied to the manual segments (-0.2 for healthy1)
                    shiftZVCSettings.toleranceGap = 0.0; % how far back should the ZVC be before we ignore it (0.5 for helathy1)
                    shiftZVCSettings.segmentProfileMod = 'half'; % should we split the window into half segments based on ZVC?
                    shiftZVCSettings.cropAllowLength = 0.5; % for data that we're only looking at half of, we're crop out this much from the segment edges, if possible. 0 to disable
                    shiftZVCSettings.removeSegmentOverflowFlag = 1; % if a manual segment exceeds the cropped data, we remove it (1) or not (0)
                    
                    % apply the half segment points policy is necessary
                    ekfDataSegTrain.time = feature_full.t; % need to be rotated for the next functions
                    ekfDataSegTrain.jointAngles = feature_full.q;
                    ekfDataSegTrain.jointVelo = feature_full.dq;
                    
                    segmentInfo_double = SegmentShiftToZVC_applyHalfSeg(segmentInfo, ekfDataSegTrain, shiftZVCSettings);
                    segmentInfo = segmentInfo_double;
                    
                    segmentInfo.timeStart = segmentInfo.startTimeVal;
                    segmentInfo.timeEnd = segmentInfo.endTimeVal;
            end
            
        case 'squats_tuat_2015'
            switch mansegMode
                case 'fullseg'
                    % combine the segments instead
                    shiftZVCSettings.offsetVal = -0.0; % initial offset applied to the manual segments (-0.2 for healthy1)
                    shiftZVCSettings.toleranceGap = 0.0; % how far back should the ZVC be before we ignore it (0.5 for helathy1)
                    shiftZVCSettings.segmentProfileMod = 'merge'; % should we split the window into half segments based on ZVC?
                    shiftZVCSettings.cropAllowLength = 0.5; % for data that we're only looking at half of, we're crop out this much from the segment edges, if possible. 0 to disable
                    shiftZVCSettings.removeSegmentOverflowFlag = 1; % if a manual segment exceeds the cropped data, we remove it (1) or not (0)
                    
                    % apply the half segment points policy is necessary
                    ekfDataSegTrain.time = feature_full.t; % need to be rotated for the next functions
                    ekfDataSegTrain.jointAngles = feature_full.q;
                    ekfDataSegTrain.jointVelo = feature_full.dq;
                    
                    segmentInfo_double = SegmentShiftToZVC_applyHalfSeg(segmentInfo, ekfDataSegTrain, shiftZVCSettings);
                    segmentInfo = segmentInfo_double;
                    
                    segmentInfo.timeStart = segmentInfo.startTimeVal;
                    segmentInfo.timeEnd = segmentInfo.endTimeVal;
            end
            
        case {'healthy1ankle', 'healthy1hip'}
            switch mansegMode
                case 'halfseg'
                    % if it's fatigue dataset, we want to create intermed segments
                    shiftZVCSettings.offsetVal = -0.0; % initial offset applied to the manual segments (-0.2 for healthy1)
                    shiftZVCSettings.toleranceGap = 0.0; % how far back should the ZVC be before we ignore it (0.5 for helathy1)
                    shiftZVCSettings.segmentProfileMod = 'half'; % should we split the window into half segments based on ZVC?
                    shiftZVCSettings.cropAllowLength = 0.5; % for data that we're only looking at half of, we're crop out this much from the segment edges, if possible. 0 to disable
                    shiftZVCSettings.removeSegmentOverflowFlag = 1; % if a manual segment exceeds the cropped data, we remove it (1) or not (0)
                    
                    % apply the half segment points policy is necessary
                    ekfDataSegTrain.time = feature_full.t; % need to be rotated for the next functions
                    ekfDataSegTrain.jointAngles = feature_full.q;
                    ekfDataSegTrain.jointVelo = feature_full.dq;
                    
                    segmentInfo_double = SegmentShiftToZVC_applyHalfSeg(segmentInfo, ekfDataSegTrain, shiftZVCSettings);
                    segmentInfo = segmentInfo_double;
                    
                    segmentInfo.timeStart = segmentInfo.startTimeVal;
                    segmentInfo.timeEnd = segmentInfo.endTimeVal;
            end
    end
    
%     % plot 
%     plot_segmentQuarter;
%     saveas(h25, fullfile(outputPath, [currInstName '_fig25_seg_' nowStr '.fig']));
%     saveas(h25, fullfile(outputPath, [currInstName '_fig25_seg_' nowStr '.png']));
%     close(h25);
    
    % converts the segment windows into 0 and 1
    segLabelAll = zeros(1, length(feature_full.t));
    for ind_windowCount = 1:1:length(segmentInfo.timeStart)
        [startVal, startInd] = findClosestValue(segmentInfo.timeStart(ind_windowCount), feature_full.t);
        [endVal, endInd] = findClosestValue(segmentInfo.timeEnd(ind_windowCount), feature_full.t);
        
        segLabelAll(startInd:endInd) = 1;
    end
    
    convertArray = [segmentInfo.timeStart; segmentInfo.timeEnd];
    for ind_windowCount = 1:1:length(convertArray)
        [closeVal, closeInd] = findClosestValue(convertArray(ind_windowCount), feature_full.t);
        startInd = closeInd - offset;
        endInd = closeInd + offset;
        
        if startInd < 1
            startInd = 1;
        end
        
        if endInd > length(feature_full.t)
            endInd = length(feature_full.t);
        end
        
        segLabelAll(startInd:endInd) = 0;
    end
    
    clusterGroupingCount = 5;
    segClusterCount = floor(length(segmentInfo.timeStart)/clusterGroupingCount);

    % not a 2011 squat fatigue file
    if length(segmentInfo.timeStart)/clusterGroupingCount ~= segClusterCount
        segClusterCount = 0;
    end
        
    segClusterToUse = [0 0];
%     segClusterToUse = [0 0; 1 length(segmentInfo.timeStart)];
%     for ind_clusterSetup = 1:segClusterCount
%         newSegClusterToAdd = [(ind_clusterSetup-1)*clusterGroupingCount+1 (ind_clusterSetup)*clusterGroupingCount];
%         segClusterToUse = [segClusterToUse; newSegClusterToAdd];
%     end
    
    for ind_clusterSetup = 1:size(segClusterToUse, 1)
        currClusterToUse = segClusterToUse(ind_clusterSetup, :);
        
%         currSegInd = currClusterToUse(1);
%         totalInd = [];
%         for ind_intermed = 1:diff(currClusterToUse)
%             [startVal, startInd] = findClosestValue(segmentInfo.timeStart(currSegInd), feature_full.t);
%             [endVal, endInd] = findClosestValue(segmentInfo.timeEnd(currSegInd), feature_full.t);
%             
%             totalInd = [totalInd startInd:endInd];
%             currSegInd = currSegInd + 1;
%         end

    % calculate the seg length and std
    metricsToMeasure{1}.name = 'seglen';
    metricsToMeasure{1}.value = mean(segmentInfo.timeEnd - segmentInfo.timeStart);
    metricsToMeasure{2}.name = 'segvar';
    metricsToMeasure{2}.value = std(segmentInfo.timeEnd - segmentInfo.timeStart);
    
    segGapCounter = [];
    for ind_seglen = 1:length(segmentInfo.timeStart) - 1
        segGapCounter(ind_seglen) = segmentInfo.timeEnd(ind_seglen+1) - segmentInfo.timeStart(ind_seglen);
    end
    
    metricsToMeasure{3}.name = 'interseglen';
    metricsToMeasure{3}.value = mean(segGapCounter);
    metricsToMeasure{4}.name = 'intersegvar';
    metricsToMeasure{4}.value = std(segGapCounter);

        if currClusterToUse(1) > 0
            [startVal, startInd] = findClosestValue(segmentInfo.timeStart(currClusterToUse(1)), feature_full.t);
            [endVal, endInd] = findClosestValue(segmentInfo.timeEnd(currClusterToUse(2)), feature_full.t);
            totalInd = startInd:endInd;
        else
            
            [startVal, startInd] = findClosestValue(segmentInfo.timeStart(1)-0.2, feature_full.t);
            [endVal, endInd]     = findClosestValue(segmentInfo.timeEnd(end)+0.2, feature_full.t);
            totalInd = startInd:endInd;
%             totalInd = 1:length(feature_full.t);
        end
        
        segLabel = segLabelAll(totalInd);

        for ind_dof = 1:size(avgWeightArray_belowThres, 2)
            weightComp = avgWeightArray_belowThres(totalInd, ind_dof);
            
            [accMax(ind_dof), thresMax(ind_dof), dirMax{ind_dof}, accMore{ind_dof}, accLess{ind_dof}] = performSeg(segLabel, weightComp, feature_full, totalInd);
            accMaxArray(ind_dof) = accMax(ind_dof).acc;
        end
        
        [maxVal, maxInd] = max(accMaxArray);

        weightCounter = zeros(size(cost_function_names));
        weightCounter(maxInd) =1 ;
        
        segLabelCountp1 = length(find(segLabel == 1));
        segLabelCountp0 = length(find(segLabel == 0));
        
        rmse_report.clusterRange = currClusterToUse;
        rmse_report.accMax = accMax;
        rmse_report.thresMax = thresMax;
        rmse_report.dirMax = dirMax;
        rmse_report.dofMax = maxInd;
        writeToOverallCSV_seg(matFilesToLoad, masterOverall, currDataset, currInstName, 0, runSettings, rmse_report, cost_function_names, param, comments, thresholdMultiplier, constThres, weightCounter, metricsToMeasure);
        writeToOverallCSV_seg(matFilesToLoad, masterSaveOverallSeg, currDataset, currInstName, 0, runSettings, rmse_report, cost_function_names, param, comments, thresholdMultiplier, constThres, weightCounter, metricsToMeasure);
        
        totalIndWeighting = length(find(totalInd > 0));
        
        %     figure; 
%     plot(ekfDataSegTrain.time, ekfDataSegTrain.jointAngles);
%     plotBoxes(gcf, segmentInfo_double.startTimeVal, segmentInfo_double.endTimeVal);
           

    if plotFig
        h = figure;

        for ind_costfunction = 1:length(cost_function_names)
            weightComp = avgWeightArray_belowThres(totalInd, ind_costfunction);
            algLabel = zeros(1, length(totalInd));
            if strcmpi(dirMax{ind_costfunction}, 'weightComp > thres')
                algLabel(weightComp > thresMax(ind_costfunction)) = 1;
            else
                algLabel(weightComp < thresMax(ind_costfunction)) = 1;
            end

            accMaxSol = accCalc(segLabel, algLabel);

            % once we've figured out the best, label the data and plot
            clf;
            subplot(211);
            plot(feature_full.t, feature_full.q);
            hold on

            ylimVal = ylim;
            minY = ylimVal(1) * 1.0005;
            maxY = ylimVal(2) * 1.0005;

            testingLabelSlide = segLabel;
            classifierLabelSlide = algLabel;

            % where to place the labels?
            labelRange = maxY - minY;
            testingLabel_max = maxY + 0.1;
            classLabel_max = maxY - 0.05*labelRange + 0.1;
            testingLabel_min = minY + 0.05*labelRange - 0.1;
            classLabel_min = minY - 0.1;

            classifierLabelSlide(classifierLabelSlide == 1) = classLabel_max;
            testingLabelSlide(testingLabelSlide == 1) = testingLabel_max;
            classifierLabelSlide(classifierLabelSlide == 0) = classLabel_min;
            testingLabelSlide(testingLabelSlide == 0) = testingLabel_min;

            plot(feature_full.t(totalInd), testingLabelSlide, 'b.');
            plot(feature_full.t(totalInd), classifierLabelSlide, 'r.');
            plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'b', 0, maxY, minY);

            ylim([classLabel_min-0.1 testingLabel_max+0.1]);


            title(['(b) gndtruth, (r) alg, segacc = ' num2str(accMaxSol.acc)]);

            subplot(212);
            plot(feature_full.t(totalInd), weightComp, 'r');
            plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'b', 0, 1.05, -0.05);
            hold on;
            plot([xlim], [thresMax(ind_costfunction) thresMax(ind_costfunction)], 'k');
            title(['Thres on ' cost_function_names{ind_costfunction} ' = ' num2str(thresMax(ind_costfunction))]);
            ylim([0 1]);

            saveas(h, fullfile(outputPath, [currInstName '_fig26_segacc_' cost_function_names{ind_costfunction} '.fig']));
            saveas(h, fullfile(outputPath, [currInstName '_fig26_segacc_' cost_function_names{ind_costfunction} '.png']));

            if ind_costfunction == maxInd
                saveas(h, fullfile(outputPath, [currInstName '_fig26_segacc_' cost_function_names{ind_costfunction} '_maxacc.fig']));
                saveas(h, fullfile(outputPath, [currInstName '_fig26_segacc_' cost_function_names{ind_costfunction} '_maxacc.png']));
                saveas(h, fullfile(masterFigureSave, [currInstName '_fig26_segacc_' cost_function_names{ind_costfunction} '_maxacc.fig']));
                saveas(h, fullfile(masterFigureSave, [currInstName '_fig26_segacc_' cost_function_names{ind_costfunction} '_maxacc.png']));
            end
        end

        close(h);

    end
    end
end

function [accMax, thresMax, dirMax, accMoreStruct, accLessStruct] = performSeg(segLabel, weightComp, feature_full, totalInd)
    testArray = 0:0.01:1;
    for i = 1:length(testArray)
        thres = testArray(i);

        algLabel = zeros(1, length(feature_full.t(totalInd)));
        algLabel(weightComp < thres) = 1;
        
        accMoreStruct(i) = accCalc(segLabel, algLabel);
        accMore(i) = accMoreStruct(i).acc;
    end
    
    for i = 1:length(testArray)
        thres = testArray(i);

        algLabel = zeros(1, length(feature_full.t(totalInd)));
        algLabel(weightComp > thres) = 1;

        accLessStruct(i) = accCalc(segLabel, algLabel);
        accLess(i) = accLessStruct(i).acc;
    end
    
    [maxMoreVal, maxMoreInd] = max(accMore);
    [maxLessVal, maxLessInd] = max(accLess);
    
    if maxMoreVal > maxLessVal
        dirMax = 'weightComp < thres';
        accMax = accMoreStruct(maxMoreInd);
        thresMax = testArray(maxMoreInd);
    else
        dirMax = 'weightComp > thres';
        accMax = accLessStruct(maxLessInd);
        thresMax = testArray(maxLessInd);
    end
end

function accMeasure = accCalc(segLabel, algLabel)
        tp = sum(and(segLabel == 1, algLabel == 1));
        fn = sum(and(segLabel == 1, algLabel == 0));
        fp = sum(and(segLabel == 0, algLabel == 1));
        tn = sum(and(segLabel == 0, algLabel == 0));
        
acc = 0.5*(tp/(tp+fn)) + 0.5*(tn/(tn+fp)); % balAcc
% acc = sum(segLabel == algLabel)/length(segLabel);

accMeasure.acc = acc;
accMeasure.tp = tp;
accMeasure.fn = fn; 
accMeasure.fp = fp;
accMeasure.tn = tn;
end