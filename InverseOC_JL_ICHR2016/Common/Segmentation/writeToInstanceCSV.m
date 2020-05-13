function writeToInstanceCSV(exportPath, subjectInfo, dataSelect, dimReductSelect, classifierSelect, aggregatorSelect, score, timeElp, observationLength, outputValues, message)
    if ~exist(exportPath, 'file')
        % if doesn't exist, write header
        header = ['runTime,dataName,dimReductName,classifierName,aggregatorName,'...
            'TrainingDatasetName,TrainingSubjectNumber,TrainingSessionNumber,TrainingExerciseName,' ...
            'TestingDatasetName,TestingSubjectNumber,TestingSessionNumber,TestingExerciseName,TestingPrimitiveName,' ...
            'TPTotal,TNTotal,TP,TN,FP,FN,TPPercentage,TNPercentage,F1_Seg,F1_Class,Acc,bAcc,MCC,AUC,' ...
            'obsLengthS,timeElp,obsLengthInd,' ...
            'normVal_q, maxInd_q1, maxInd_q2, normVal_q_all, maxInd_q1_all, maxInd_q2_all,' ...
            'normVal_dq, maxInd_dq, maxInd_dq2, normVal_dq_all, maxInd_dq1_all, maxInd_dq2_all'];
    else
        header = '';
    end
    
    if ~exist('outputValues', 'var') || isempty(outputValues)
        outputValues.totalLength = 0;
        outputValues.normVal_q = 0;
        outputValues.maxInd_q1 = 0;
        outputValues.maxInd_q2 = 0;
        outputValues.normVal_q_all = 0;
        outputValues.maxInd_q1_all = 0;
        outputValues.maxInd_q2_all = 0;
        outputValues.normVal_dq = 0;
        outputValues.maxInd_dq1 = 0;
        outputValues.maxInd_dq2 = 0;
        outputValues.normVal_dq_all = 0;
        outputValues.maxInd_dq1_all = 0;
        outputValues.maxInd_dq2_all = 0;
    end
    
    
    if ~exist('message', 'var')
        message = ''; 
    end
    
    fId = fopenCheck(exportPath, header);
    
    if exist('subjectInfo', 'var')
        % parse the subjectInfo
        if isfield(subjectInfo, 'testing_datasetName')
            
            % build the px and exercise name
            if ischar(subjectInfo.training_subjectNumber)
                pxString = subjectInfo.training_subjectNumber;
                sessString = subjectInfo.training_sessionNumber;
                exerciseString = subjectInfo.training_exerciseName;
            else
                % build the px name
                pxString = [];
                for i = 1:length(subjectInfo.training_subjectNumber)
                    pxString = [pxString '-' num2str(subjectInfo.training_subjectNumber(i))];
                end
                pxString = pxString(2:end); % remove the first dash
                
                % build the sess name
                sessString = [];
                for i = 1:length(subjectInfo.training_sessionNumber)
                    sessString = [sessString '-' num2str(subjectInfo.training_sessionNumber(i))];
                end
                sessString = sessString(2:end); % remove the first dash
                
                 % build the exercise
                exerciseString = [];
                if iscell(subjectInfo.training_exerciseName)
                    trainingExerciseNameSorted = sort(subjectInfo.training_exerciseName);
                    for i = 1:length(subjectInfo.training_exerciseName)
                        exerciseString = [exerciseString '-' trainingExerciseNameSorted{i}];
                    end
                    exerciseString = exerciseString(2:end); % remove the first dash
                else
                    exerciseString = subjectInfo.training_exerciseName;
                end
            end

            trainingDatasetName = subjectInfo.training_datasetName;
            trainingSubjectNumber = pxString;
            trainingSessionNumber = sessString;
            trainingExerciseName = exerciseString;
            
            testingDatasetName = subjectInfo.testing_datasetName;
            testingSubjectName = num2str(subjectInfo.testing_subjectNumber); % the testing should always only consist of a single px/sess
            testingSessionNumber = num2str(subjectInfo.testing_sessionNumber);
            testingExerciseName = subjectInfo.testing_exerciseName;
            testingPrimitiveName = subjectInfo.testing_primitiveName;
        else % legacy support
            trainingDatasetName = '-';
            trainingSubjectNumber = '-';
            trainingSessionNumber = '-';
            trainingExerciseName = '-';
            
            testingDatasetName = subjectInfo.datasetName;
            testingSubjectName = num2str(subjectInfo.subjectNumber);
            testingSessionNumber = num2str(subjectInfo.sessionNumber);
            testingExerciseName = subjectInfo.exerciseName;
            testingPrimitiveName = '-';
        end
        
        score = checkScore(score);
        
        fprintf(fId, '%s,%s,%s,%s,%s ,%s,%s,%s,%s, %s,%s,%s,%s,%s ,%u,%u,%u,%u ,%u,%u,%1.3f,%1.3f, %1.3f,%1.3f,%1.3f,%1.3f,%1.3f,%1.3f ,%f,%f,%u,  %1.3f,%u,%u,%1.3f,%u,%u, %1.3f,%u,%u,%1.3f,%u,%u ,%s \n', ...
            datestr(now, 'YYYY-mm-DD hh-MM-ss'), ...
            dataSelect.settings.settingName, dimReductSelect.settingName, classifierSelect.settingName, aggregatorSelect.settingName, ...
            trainingDatasetName, trainingSubjectNumber, trainingSessionNumber, trainingExerciseName, ...
            testingDatasetName, testingSubjectName, testingSessionNumber, testingExerciseName, testingPrimitiveName, ...
            score.segmentTruePosTotal, score.segmentTrueNegTotal, score.segmentTruePos, score.segmentTrueNeg, ...
            score.segmentFalsePos, score.segmentFalseNeg, score.segmentTruePosPercentage, score.segmentTrueNegPercentage, ...
            score.fScore_Seg, score.fScore_Class, score.accuracy, score.bAcc, score.MCC, score.AUC, ...
            observationLength, timeElp, outputValues.totalLength, ...
            outputValues.normVal_q, outputValues.maxInd_q1, outputValues.maxInd_q2, ...
            outputValues.normVal_q_all, outputValues.maxInd_q1_all, outputValues.maxInd_q2_all, ...
            outputValues.normVal_dq, outputValues.maxInd_dq1, outputValues.maxInd_dq2, ...
            outputValues.normVal_dq_all, outputValues.maxInd_dq1_all, outputValues.maxInd_dq2_all, ...
            message);
    end
    
    fclose(fId);
end

function score = checkScore(score)
    % check to make sure all the fields are populated. If not, fill them in
    % with zero
    
    if ~isfield(score, 'segmentTruePosTotal')
        score.segmentTruePosTotal = 0;
    end
    
    if ~isfield(score, 'segmentTrueNegTotal')
        score.segmentTrueNegTotal = 0;
    end
    
    if ~isfield(score, 'segmentTruePos')
        score.segmentTruePos = 0;
    end
    
    if ~isfield(score, 'segmentTrueNeg')
        score.segmentTrueNeg = 0;
    end
    
    if ~isfield(score, 'segmentFalsePos')
        score.segmentFalsePos = 0;
    end
    
    if ~isfield(score, 'segmentFalseNeg')
        score.segmentFalseNeg = 0;
    end
    
    if ~isfield(score, 'segmentTruePosPercentage')
        score.segmentTruePosPercentage = 0;
    end
    
    if ~isfield(score, 'segmentTrueNegPercentage')
        score.segmentTrueNegPercentage = 0;
    end
    
    if ~isfield(score, 'fScore_Seg')
        score.fScore_Seg = 0;
    end
    
    if ~isfield(score, 'fScore_Class')
        score.fScore_Class = 0;
    end
    
    if ~isfield(score, 'accuracy')
        score.accuracy = 0;
    end
    
    if ~isfield(score, 'bAcc')
        score.bAcc = 0;
    end
    
    if ~isfield(score, 'MCC')
        score.MCC = 0;
    end
    
    if ~isfield(score, 'AUC')
        score.AUC = 0;
    end
end