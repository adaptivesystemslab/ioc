function runSharcnet(mFileName, outSource)
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
    
    [structTrainingExerciseList, structTestingExerciseList] = exerciseSet(movementType, testingType);
    
    batchSettings.exportPathSuffix = structExportPath;
    
    batchSettings.trainingPatientList = trainingPxList;
    batchSettings.trainingExerciseList = structTrainingExerciseList;
    
    batchSettings.testingPatientList = testingPxList;
    batchSettings.testingExerciseList = structTestingExerciseList;
    
    batchSettings.dataSelectList = dataSelect;
    
    setBatchSettings;
    
    if length(fileNameParse) > 5 && ~isempty(fileNameParse{6})
        % we're passing in specific params
        batchSettings.dimReductionList = fileNameParse(6);
        
        if strcmpi(batchSettings.dimReductionList, 'ksPCA2')
            batchSettings.dimReductionList = {'ksPCA2_60_32_none_delta_1000'};
        end
    end
    
    if length(fileNameParse) > 6 && ~isempty(fileNameParse{7})
        batchSettings.aggregatorList = fileNameParse(7);
    end
    
    if length(fileNameParse) > 7 && ~isempty(fileNameParse{8})
%         batchSettings.classifierList = fileNameParse(8);
        
        if strcmpi(fileNameParse{8}, 'RBF')
            newClassifierArray{1} = 'RBF_10';
            newClassifierArray{2} = 'RBF_20';
        elseif strcmpi(fileNameParse{8}, 'kNN')
            newClassifierArray{1} = 'kNN_3';
            newClassifierArray{2} = 'kNN_9';
        elseif strcmpi(fileNameParse{8}, 'SVM')
            newClassifierArray{1} = 'SVM_radial';
            newClassifierArray{2} = 'SVM_polynomial';
            newClassifierArray{3} = 'SVM_linear';
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
    
    
%     batchSettings.batchInstancePath = ['script20_PCWINrun2_' datestr(now, 'yyyy-mm-dd')]
%     batchSettings.aggregatorList = {'None'}
%     batchSettings.classifierList = batchSettings.classifierList(1)
    
    batchSettings    
    batchClassifier_healthy1(batchSettings);

    fprintf('Completed %s\n', datestr(now));
end

function pxList = patientList(patientType)
% templatesToUse{1} = [2 3 4 6 8];
% templatesToUse{2} = [9 15 16 19 20];
% 
% templatesToUse{3} = [5 7 10 11 12]; % KHEF
% templatesToUse{4} = [1 13 14 17 18];

    switch patientType
        case 1
            pxList{1} = [2 3 4 6 8];
            
        case 2
            pxList{1} = [9 15 16 19 20];
            
        case 3
            pxList{1} = [5 7 10 11 12]; % KHEF

        case 4
            pxList{1} = [1 13 14 17 18];
            
        case 9
            for i = 1:20
                pxList{i} = i;
            end
            
        case 11
            pxList{1} = [2 3 4 6 8];
            pxList{2} = [9 15 16 19 20];
            pxList{3} = [5 7 10 11 12]; % KHEF
            pxList{4} = [1 13 14 17 18];
            
        case 12
            a = 1:20;
            pxList{1} = setxor(a, [2 3 4 6 8]);
            pxList{2} = setxor(a, [9 15 16 19 20]);
            pxList{3} = setxor(a, [5 7 10 11 12]); % KHEF
            pxList{4} = setxor(a, [1 13 14 17 18]);
            
        case 13
            pxList{2} = [2 3 4 6 8];
            pxList{1} = [9 15 16 19 20];
            pxList{4} = [5 7 10 11 12]; % KHEF
            pxList{3} = [1 13 14 17 18];
            
            
            
        case 21
            pxList{1} = [2 3 4 6 8];
            pxList{2} = [9 15 16 19 20];
            
        case 22
            pxList{2} = [2 3 4 6 8];
            pxList{1} = [9 15 16 19 20];
            
        case 31
            pxList{1} = [5 7 10 11 12]; % KHEF
            pxList{2} = [1 13 14 17 18];
            
        case 32
            pxList{2} = [5 7 10 11 12]; % KHEF
            pxList{1} = [1 13 14 17 18];
            
        otherwise
            % use negative numbers to denote the actual px
            pxList{1} = -1 * patientType;
    end
end

function [structTrainingExerciseList, structTestingExerciseList] = exerciseSet(movementType, testingType)
    
    if strcmpi(movementType, 'All')
        structTrainingExerciseList = {'HFEO_SUP_SLO1', 'KEFO_SIT_SLO1', 'KHEF_SUP_SLO1', 'SQUA_STD_SLO1', 'STSO_SIT_SLO1', ...
            'HFEO_SUP_SLO2', 'KEFO_SIT_SLO2', 'KHEF_SUP_SLO2', 'SQUA_STD_SLO2', 'STSO_SIT_SLO2'};
        structTestingExerciseList = {'HFEO_SUP_SLO1', 'KEFO_SIT_SLO1', 'KHEF_SUP_SLO1', 'SQUA_STD_SLO1', 'STSO_SIT_SLO1', ...
            'HFEO_SUP_SLO2', 'KEFO_SIT_SLO2', 'KHEF_SUP_SLO2', 'SQUA_STD_SLO2', 'STSO_SIT_SLO2'};
        
    elseif strcmpi(movementType(1:5), 'Test2')
        structTrainingExerciseList = {'HFEO_SUP_SLO1', 'KEFO_SIT_SLO1', 'KHEF_SUP_SLO1', 'SQUA_STD_SLO1', 'STSO_SIT_SLO1'};
        structTestingExerciseList = {'HFEO_SUP_SLO2', 'KEFO_SIT_SLO2', 'KHEF_SUP_SLO2', 'SQUA_STD_SLO2', 'STSO_SIT_SLO2'};

    elseif strcmpi(movementType(1:5), 'Test3')
        
        structTestingExerciseList = {'HFEO_SUP_SLO2', 'KEFO_SIT_SLO2', 'KHEF_SUP_SLO2', 'SQUA_STD_SLO2', 'STSO_SIT_SLO2'};

        
        switch str2num(movementType(7:8))
            case 11
                structTrainingExerciseList = {'KEFO_SIT_SLO1', 'STSO_SIT_SLO1'};
            case 12
                structTrainingExerciseList = {'SQUA_STD_SLO1'};
            case 13
                structTrainingExerciseList = {'HFEO_SUP_SLO1', 'KHEF_SUP_SLO1'};
            case 21
                structTrainingExerciseList = {'KEFO_SIT_SLO1', 'HFEO_SUP_SLO1', 'KHEF_SUP_SLO1'};
            case 22
                structTrainingExerciseList = {'STSO_SIT_SLO1', 'SQUA_STD_SLO1'};
            case 31
                structTrainingExerciseList = {'KEFO_SIT_SLO1'};
            case 32
                structTrainingExerciseList = {'STSO_SIT_SLO1', 'KHEF_SUP_SLO1', 'SQUA_STD_SLO1'};
            case 33
                structTrainingExerciseList = {'HFEO_SUP_SLO1'};
            case 40
                structTrainingExerciseList = {'HFEO_SUP_SLO1', 'KEFO_SIT_SLO1', 'KHEF_SUP_SLO1', 'SQUA_STD_SLO1', 'STSO_SIT_SLO1'};
        end
    else
        if isempty(str2num(movementType(end)))
            structTrainingExerciseList{1} = [movementType num2str(1)];
            structTestingExerciseList{1} = [movementType num2str(2)];
        else
            structTrainingExerciseList{1} = [movementType];
            structTestingExerciseList{1} = [movementType];
        end
        
        if strcmpi(testingType, 'Test1')
            % use everything for the testing sample
            structTestingExerciseList = {'HFEO_SUP_SLO2', 'KEFO_SIT_SLO2', 'KHEF_SUP_SLO2', 'SQUA_STD_SLO2', 'STSO_SIT_SLO2'};
        end
        
%         structTestingExerciseList = structTrainingExerciseList;
    end
end