function runSharcnet_TRO(mFileName, outSource)
    fprintf('-------------------------------------------------------------\n');
    fprintf('-------------------------------------------------------------\n');
    fprintf('Starting SHARCNET-Segmentation\n');
    fprintf('Started %s\n', datestr(now));
    
    batchSettings.mFileName = mFileName;
    batchSettings.batchInstancePath = outSource;
    
    fileNameParse = regexp(mFileName,'z','split');
%     fileNameParse = strsplit('z', mFileName);
    
    movementType = fileNameParse{2};
    dataSelect = fileNameParse(3);
    trainingPatientType = str2num(fileNameParse{4});
    testingPatientType = str2num(fileNameParse{5});
    
    if length(fileNameParse) > 8
    testingType = fileNameParse{9};
    else
        testingType = '';
    end
    
    % generate some lists
%     structExportPath = [datestr(now, 'YYYY-mm-DD-hh-MM-ss') '-' movementType '-' dataSelect{1}];
%     structExportPath = [movementType '-' dataSelect{1}];
    structExportPath = [dataSelect{1}];
    trainingPxList = patientList(trainingPatientType);
    testingPxList = patientList(testingPatientType);
    
    [structTrainingExerciseList, structTestingExerciseList, exerciseListFull, fullSegmentsFlag] = exerciseSet(movementType);
    
    batchSettings.exportPathSuffix = structExportPath;
    
    batchSettings.trainingPatientList = trainingPxList;
    batchSettings.trainingExerciseList = structTrainingExerciseList;
    
    batchSettings.testingPatientList = testingPxList;
    batchSettings.testingExerciseList = structTestingExerciseList;
    batchSettings.exerciseListFull = exerciseListFull;
    batchSettings.fullSegmentsFlag = fullSegmentsFlag;
    
    batchSettings.dataSelectList = dataSelect;
    
    setBatchSettings;
    
    if length(fileNameParse) > 5 && ~isempty(fileNameParse{6})
        % we're passing in specific params
        batchSettings.dimReductionList = fileNameParse(6);
    end
    
    if length(fileNameParse) > 6 && ~isempty(fileNameParse{7})
        batchSettings.aggregatorList = fileNameParse(7);
    end
    
    if length(fileNameParse) > 7 && ~isempty(fileNameParse{8})
%         batchSettings.classifierList = fileNameParse(8);
        
        if strcmpi(fileNameParse{8}, 'RBF')
            newClassifierArray{1} = 'RBF_10';
            newClassifierArray{2} = 'RBF_20';
        elseif strcmpi(fileNameParse{8}, 'SVM')
            newClassifierArray{1} = 'SVM_linear';
            newClassifierArray{2} = 'SVM_polynomial';
            newClassifierArray{3} = 'SVM_radial';
        elseif strcmpi(fileNameParse{8}, 'SVM_poly')
            newClassifierArray{1} = 'SVM_polynomial';
        elseif strcmpi(fileNameParse{8}, 'NN')
            newClassifierArray{1} = 'NN_10_10';
            newClassifierArray{2} = 'NN_10_10_10';
            newClassifierArray{3} = 'NN_20_20_20'; % 'NN_10_10', 'NN_10_10_10', 'NN_20_10_5', 'NN_20_20_20'
        else
            newClassifierArray = fileNameParse(8);
        end
        
          batchSettings.classifierList = newClassifierArray;
          
%             classifierStack = {'kNN_3', 'kNN_9', 'QDA', ...
%         'RBF_10', 'RBF_20', ...
%         'SVM_linear', 'SVM_polynomial', 'SVM_radial'...
%         'NN_10_10'};
    end
    
    
    batchSettings    
    batchClassifier_TRO(batchSettings);

    fprintf('Completed %s\n', datestr(now));
end

function pxList = patientList(patientType)
% templatesToUse{1} = [2 3 4 6 8];
% templatesToUse{2} = [9 15 16 19 20];
% 
% templatesToUse{3} = [5 7 10 11 12]; % KHEF
% templatesToUse{4} = [1 13 14 17 18];

    switch patientType            
        otherwise
            pxList{1} = [1];
    end
end

function [structTrainingExerciseList, structTestingExerciseList, exerciseListFull, fullSegmentsFlag] = exerciseSet(movementType)

    movementTypeSplit = strsplit('_', movementType);
    
    
    switch movementTypeSplit{1}
        case 'Test3'
            
            exerciseListFull = {'BAD', 'BAL180', 'BAL90', 'BAR180', 'BAR90', 'BAU', 'LAL180', ...
                'LAL90', 'LAR180', 'LAR90', 'LER', 'LKAL', 'LKAR', 'LKE', 'LKR', 'LPAL', 'LPAR', ...
                'LPE', 'LPR', 'LPUAD', 'LPUE', 'LPUR', 'MLL', 'MLR', 'MRL', 'MRR', 'RAL180', ...
                'RAL90', 'RAR180', 'RAR90', 'RKAL', 'RKAR', 'RKE', 'RKR', 'RPAL', 'RPAR', 'RPE', ...
                'RPR', 'RRL', 'SQD', 'SQU', 'WLL', 'WLR', 'WRL', 'WRR'};
            
            fullSegmentsFlag = 0;
            structTestingExerciseList = {'ALLMOTIONS'};
            
            switch str2num(movementTypeSplit{2})
                case 11 % right arm
                    structTrainingExerciseList = {'RAL180', 'RAL90', 'RAR180', 'RAR90', 'RPE', 'RPR'};
                case 12 % left arm
                    structTrainingExerciseList = {'LAL180', 'LAL90', 'LRAR180', 'LRAR90', 'LRPE', 'LRPR'};
                case 13 % right leg
                    structTrainingExerciseList = {'LKE', 'LKR'};
                case 14 % left leg
                    structTrainingExerciseList = {'RKE', 'RKR'};
                case 15 % upper body
                    structTrainingExerciseList = {'BAD', 'BAL180', 'BAL90', 'BAR180', 'BAR90', 'BAU', ...
                        'LKAL', 'LKAR', 'LPAL', 'LPAR', 'LPUAD', 'LPUE', 'LPUE', 'LPUE', ...
                        'LPUR', 'RKAL', 'RKAR', 'RPAL', 'RPAR'};
                    %             case 16 % lower body
                    %                 structTrainingExerciseList = {'', 'SQUA_STD_SLO1'};
                case 17 % full body
                    structTrainingExerciseList = {'MLL', 'MLR', 'MRL', 'MRR', 'WLL', 'WLR', 'WRL', 'WRR'};
                    
                case 21
                    structTrainingExerciseList = {'BAD', 'BAU'};
                case 22
                    structTrainingExerciseList = {'BAL180', 'BAR180'};
                case 23
                    structTrainingExerciseList = {'BAL90', 'BAR90'};
                case 24
                    structTrainingExerciseList = {'LAL180', 'LAR180'};
                case 25
                    structTrainingExerciseList = {'LAL90', 'LAL90'};
                case 26
                    structTrainingExerciseList = {'LKAR', 'LKAR', 'LKE', 'LKR'};
                case 27
                    structTrainingExerciseList = {'LPAL', 'LPAR', 'LPE', 'LPR', 'LPUAD', 'LPUE', 'LPUR'};
                case 28
                    structTrainingExerciseList = {'MLL', 'MLR', 'MRL', 'MRR'};
                case 29
                    structTrainingExerciseList = {'RAL180', 'RAR180'};
                case 30
                    structTrainingExerciseList = {'RAL90', 'RAR90'};
                case 31
                    structTrainingExerciseList = {'RKAL', 'RKAR', 'RKE', 'RKR'};
                case 32
                    structTrainingExerciseList = {'RPAL', 'RPAR', 'RPE', 'RPR'};
                case 33
                    structTrainingExerciseList = {'SQU', 'SQD'};
                case 34
                    structTrainingExerciseList = {'WLL', 'WLR', 'WRL', 'WRR'};
                    
                    
            end
            
            %         case 'Reduced_All'
            %             structTrainingExerciseList = {'RAR90', 'LPE', 'LKE', 'LAR90', 'BAR90', 'BAR180'};
            %             structTestingExerciseList = {'ALLMOTIONS'};
            
            %         case 'Larger_Limited'
            %             structTrainingExerciseList = {'RAR90', 'LPE', 'LKE', 'LAR90', 'BAR90', 'BAR180', 'BAD', 'MLR', 'SQD'};
            %             structTestingExerciseList = {'RAR90', 'LPE', 'LKE', 'LAR90', 'BAR90', 'BAR180', 'BAD', 'MLR', 'SQD'};
            
            %         case 'Larger_All'
            % %             structTrainingExerciseList = {'RAR90', 'LPE', 'LKE', 'LAR90', 'BAR90', 'BAR180', 'BAD', 'MLR', 'SQD'};
            %             structTrainingExerciseList = {'RAR90', 'RAL90', 'LPE', 'LPR', 'LKE', 'LKR', ...
            %                 'LAR90', 'LAL90', 'BAR90', 'BAL90', 'BAR180', 'BAL180', 'BAD', 'BAU', ...
            %                 'MLR', 'MLL', 'MRR', 'MRL', 'SQD', 'SQU'};
            %             structTestingExerciseList = {'ALLMOTIONS'};
            
        case 'HalfSeg'
            structTrainingExerciseList = {'RAR90', 'RAL90', 'LKAL' ,'LPE', 'LPR', 'LKE', 'LKR', ...
                'LAR90', 'LAL90', 'BAR90', 'BAL90', 'BAR180', 'BAL180', 'BAD', 'BAU', ...
                'MLR', 'MLL', 'MRR', 'MRL', 'SQD', 'SQU'};
            
            exerciseListFull = {'BAD', 'BAL180', 'BAL90', 'BAR180', 'BAR90', 'BAU', 'LAL180', ...
                'LAL90', 'LAR180', 'LAR90', 'LER', 'LKAL', 'LKAR', 'LKE', 'LKR', 'LPAL', 'LPAR', ...
                'LPE', 'LPR', 'LPUAD', 'LPUE', 'LPUR', 'MLL', 'MLR', 'MRL', 'MRR', 'RAL180', ...
                'RAL90', 'RAR180', 'RAR90', 'RKAL', 'RKAR', 'RKE', 'RKR', 'RPAL', 'RPAR', 'RPE', ...
                'RPR', 'RRL', 'SQD', 'SQU', 'WLL', 'WLR', 'WRL', 'WRR'};
            
            fullSegmentsFlag = 0;
            structTestingExerciseList = {'ALLMOTIONS'};
            
        case 'FullSeg'
            structTrainingExerciseList = {'RAR90', 'LPE', 'LKE', 'LAR90', 'BAR90', 'BAR180', 'BAD', 'MLR', 'SQD'};
            
            exerciseListFull = {'BAD', 'BAR180', 'BAR90', 'LAR180', 'LAR90', 'LKAR', ...
                'LKE', 'LPAR', 'LPE', 'LPUE', 'MLR', 'RAR180', 'RAR90', 'RKAR', 'RKE', 'RPAR', 'RPE', 'SQD', 'WLR'}; % full list
            
            fullSegmentsFlag = 1;
            structTestingExerciseList = {'ALLMOTIONS'};

        otherwise
            structTrainingExerciseList{1} = [movementType];
            
            exerciseListFull = {'BAD', 'BAL180', 'BAL90', 'BAR180', 'BAR90', 'BAU', 'LAL180', ...
                'LAL90', 'LAR180', 'LAR90', 'LER', 'LKAL', 'LKAR', 'LKE', 'LKR', 'LPAL', 'LPAR', ...
                'LPE', 'LPR', 'LPUAD', 'LPUE', 'LPUR', 'MLL', 'MLR', 'MRL', 'MRR', 'RAL180', ...
                'RAL90', 'RAR180', 'RAR90', 'RKAL', 'RKAR', 'RKE', 'RKR', 'RPAL', 'RPAR', 'RPE', ...
                'RPR', 'RRL', 'SQD', 'SQU', 'WLL', 'WLR', 'WRL', 'WRR'};
            
            fullSegmentsFlag = 0;
            structTestingExerciseList = {'ALLMOTIONS'};
    end
    
    if length(movementTypeSplit) > 1 && strcmpi(movementTypeSplit{2}, 'rnd')
        % randomizing the input
        randExerciseList = randperm(length(exerciseListFull), floor(length(exerciseListFull)/2));
        structTrainingExerciseList = exerciseListFull(randExerciseList);
    end
end