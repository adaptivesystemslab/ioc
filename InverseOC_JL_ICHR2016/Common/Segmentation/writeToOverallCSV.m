function writeToOverallCSV(exportPath, dataSelect, dimReductSelect, classifierSelect, aggregatorSelect, score, timeElp, obsLength, message,extraHeader,extraEntry)

    if ~exist('message', 'var')
        message = ' ';
    end
    
    if ~exist('extraHeader', 'var')
        extraHeader = '';
        extraEntry = '';
    end

    if ~exist(exportPath, 'file')
        % if doesn't exist, write header
        header = ['runTime,dataName,dimReductName,classifierName,aggregatorName,' ...
            'datasetTraining,subjsetTraining,exercisesetTraining,mansegTraining,' ...
            'datasetTesting,subjsetTesting,exercisesetTesting,mansegTesting,' ...
            'TPTotal,TNTotal,TP,TN,FP,FN,TPPercentage,TNPercentage,F1_Seg,F1_Seg_std,F1_Class,F1_Class_std,Acc,Acc_std,bAcc,bAcc_std,MCC,MCC_std,AUC,AUC_std,' ...
            'obsLength,timeTaken,confMtx (GndTruth\\Class),' extraHeader];
    else
        header = '';
    end
    
    fId = fopenCheck(exportPath, header);
    
    % build the px and exercise name 
% % %     if isfield(dataSelect, 'dataPackage') && isfield(dataSelect, 'trainingDataPackage')
% % %         [trainingDatasetString, trainingPatientString, trainingExerciseString] = ...
% % %             compileString(dataSelect.dataPackageTraining, dataSelect.exerciseStruct);
% % %         
% % %         [testingDatasetString, testingPatientString, testingExerciseString] = ...
% % %             compileString(dataSelect.dataPackage, dataSelect.exerciseStruct);
% % %     else
% % %         trainingDatasetString = '-';
% % %         trainingPatientString = '-';
% % %         trainingExerciseString = '-';
% % %         
% % %         testingDatasetString = [];
% % %         testingPatientString = [];
% % %         testingExerciseString = [];
% % %     
% % %         for i = 1:length(dataSelect.patientSelection.pxCount)
% % %             testingPatientString = [testingPatientString '-' num2str(dataSelect.patientSelection.pxCount(i))];
% % %         end
% % %         testingPatientString = testingPatientString(2:end); % remove the first dash
% % %         
% % %         
% % %         for i = 1:length(dataSelect.patientSelection.exerciseName)
% % %             testingExerciseString = [testingExerciseString '-' dataSelect.patientSelection.exerciseName{i}];
% % %         end
% % %         testingExerciseString = testingExerciseString(2:end); % remove the first dash
% % %     end

    if isfield(dataSelect, 'dataPackageTraining')
        [trainingDatasetString, trainingPatientString, trainingExerciseString, trainingManSegString] = ...
            compileString(dataSelect.dataPackageTraining, dataSelect.exerciseStructTraining);
    else
        trainingDatasetString = '-';
        trainingPatientString = '-';
        trainingExerciseString = '-';
        trainingManSegString = '-';
    end

    if isfield(dataSelect, 'testingExerciseName')
        [testingDatasetString, testingPatientString, testingExerciseString, testingManSegString] = ...
            compileString(dataSelect.dataPackage, dataSelect.exerciseStruct);
    elseif isfield(dataSelect, 'dataPackage')
        [testingDatasetString, testingPatientString, testingExerciseString, testingManSegString] = ...
            compileString(dataSelect.dataPackage, dataSelect.exerciseStruct);
    else
        testingDatasetString = [];
        testingPatientString = [];
        testingExerciseString = [];

        for i = 1:length(dataSelect.patientSelection.pxCount)
            testingPatientString = [testingPatientString '-' num2str(dataSelect.patientSelection.pxCount(i))];
        end
        testingPatientString = testingPatientString(2:end); % remove the first dash


        for i = 1:length(dataSelect.patientSelection.exerciseName)
            testingExerciseString = [testingExerciseString '-' dataSelect.patientSelection.exerciseName{i}];
        end
        testingExerciseString = testingExerciseString(2:end); % remove the first dash
        
        testingManSegString = '-';
    end

    if ~isempty(score)
        score = checkScore(score);
        fprintf(fId, '%s,%s,%s,%s,%s ,%s,%s,%s,%s, %s,%s,%s,%s, %u,%u,%u,%u, %u,%u, %1.3f,%1.3f, %1.3f,%1.3f,%1.3f,%1.3f, %1.3f,%1.3f,%1.3f,%1.3f, %1.3f,%1.3f,%1.3f,%1.3f, %3.2f,%3.2f, %s,%s \n', ...
            datestr(now, 'YYYY-mm-DD hh-MM-ss'), dataSelect.settings.settingName, dimReductSelect.settingName, classifierSelect.settingName, aggregatorSelect.settingName, ...
            trainingDatasetString, trainingPatientString, trainingExerciseString, trainingManSegString,...
            testingDatasetString, testingPatientString, testingExerciseString, testingManSegString, ...
            score.segmentTruePosTotal, score.segmentTrueNegTotal, score.segmentTruePos, score.segmentTrueNeg, ...
            score.segmentFalsePos, score.segmentFalseNeg, ...
            score.segmentTruePosPercentage, score.segmentTrueNegPercentage, ...
            score.fScore_Seg, score.fScore_Seg_std, score.fScore_Class, score.fScore_Class_std, ...
            score.accuracy, score.accuracy_std, score.bAcc, score.bAcc_std, ...
            score.MCC, score.MCC_std, score.AUC, 0, ...
            obsLength, timeElp, ...
            message, extraEntry);
        
%         for i = 1:size(score.confMtx, 1)
%             for j = 1:size(score.confMtx, 2)
%                 fprintf(fId, ',%u', score.confMtx(i, j));
%             end
%         end

%         fprintf(fId, ',%u,%u,%u,%u,%u,%u,%u,%u,%s', ...
%             score.breakdownTPSum_higProb, score.breakdownTPSum_lowProb, ...
%             score.breakdownTNSum_higProb, score.breakdownTNSum_lowProb, ...
%             score.breakdownFNSum_higProb, score.breakdownFNSum_lowProb, ...
%             score.breakdownFPSum_higProb, score.breakdownFPSum_lowProb, ...
%             message);
    
    else
        fprintf(fId, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,%s \n', ...
            datestr(now, 'YYYY-mm-DD hh-MM-ss'), ...
            dataSelect.settings.settingName, dimReductSelect.settingName, ...
            classifierSelect.settingName, aggregatorSelect.settingName, ...
            trainingDatasetString, trainingPatientString, trainingExerciseString, trainingManSegString, ...
            testingDatasetString, testingPatientString, testingExerciseString, testingManSegString, ...
            message);
    end
    
    fclose(fId);
end

function [datasetStr, patientStr, exerciseStr, mansegStr] = compileString(dataPackage, exerciseStruct)
    datasetStr = [];
    patientStr = [];
    exerciseStr = [];
    mansegStr = [];
    
    % sort all the packages first so it's in order
    localExerciseType = sort(exerciseStruct.exerciseType);

    for ind = 1:length(dataPackage)
        currDataPackage = dataPackage{ind};

        datasetStr = [datasetStr '-' currDataPackage.dataset];

        for i = 1:length(currDataPackage.patient)
            patientStr = [patientStr '-' num2str(currDataPackage.patient(i))];
        end
        
        if ~isfield(currDataPackage, 'manSeg')
            mansegStr = '';
        else
            mansegStr = [mansegStr '-' currDataPackage.manSeg];
        end
    end

    for i = 1:length(localExerciseType)
        exerciseStr = [exerciseStr '-' localExerciseType{i}];
    end
    
    datasetStr = datasetStr(2:end);
    patientStr = patientStr(2:end); % remove the first dash
    exerciseStr = exerciseStr(2:end); % remove the first dash
    mansegStr = mansegStr(2:end);
end

function score = checkScore(score)
    % check to make sure all the fields are populated. If not, fill them in
    % with zero
    
    if ~isfield(score, 'segmentTrueNegTotal')
        score.segmentTrueNegTotal = 0;
    end
    
    if ~isfield(score, 'segmentTrueNeg')
        score.segmentTrueNeg = 0;
    end
    
    if ~isfield(score, 'segmentTrueNegPercentage')
        score.segmentTrueNegPercentage = 0;
    end
end