% this script loads the content of databaseSpecTablebasespec, runs EKF on the mocap databaseSpecTable, 
% then saves the output into a target folder
clearvars
clc

runModelParam = 0;
runEkf_mocap = 1;
runAnalysis = 0;

overwriteExisting = 1;
visualize_pose = 1;
visualize_main = 1;

basepathSource = 'D:\aslab\data\FullBody_IIT_2017'; 
basepathTarget = 'D:\aslab\data_IK\FullBody_IIT_2017';

filepathModelXml = fullfile('.', 'model', 'iit_v10.xml');

% filepaths to the source databaseSpecTable and target databaseSpecTable
timeStamp = '2018_08_20_checktest';
% timeStamp = datestr(now, 'yyyy_mm_dd_HH_MM_SS');

filepathDecoder = fullfile(basepathSource,'databaseSpec.csv'); 
filepathBlacklistDecoder = fullfile(basepathSource,'modelParamBlackList.csv'); 
filepathSource = fullfile(basepathSource,filesep);
filepathTarget = fullfile(basepathTarget,timeStamp,filesep);

% setting the path for common files
addpath(genpath('.\support'));
addpath(genpath('..\common'));
addpath(genpath('..\toolboxes'));
addpath(genpath('..\ekf'));

% loading databaseSpec
databaseSpecTable = readtable(filepathDecoder);
exerciseNames = fieldnames(databaseSpecTable);
exerciseNames = exerciseNames(2:12); % remove headers and other non-exercise fields

externalParam_mocap = [];
externalParam_mocap.modalityType = 'mocap'; % mocap imu
externalParam_mocap.modalityCalibration = 'mocap'; % setting to set calibration matrices
externalParam_mocap.linkDefinition = 'X00'; % setting to set kinematic for trc model (initial X00 anth)
externalParam_mocap.modelBase = 'floating'; % floating hip
externalParam_mocap.subSuffix = '';
externalParam_mocap.externalTrcData = [];
externalParam_mocap.externalImuData = [];
externalParam_mocap.externalBaseFrameData = [];
externalParam_mocap.overwriteExisting = overwriteExisting;
externalParam_mocap.visualize_pose = visualize_pose;
externalParam_mocap.visualize_main = visualize_main;

% compile file list for file loading
indFileEntry = 0;
for ind_subj = 1:length(databaseSpecTable.SUBJECT)
    currSubjId = databaseSpecTable.SUBJECT(ind_subj);
    for ind_exercise = 1:length(exerciseNames)
        % load a specific subject and exercise name
        currExerId = exerciseNames{ind_exercise};
        currSessId = 1;
        
        % make sure that it's "subject01" instead of "subject1" 
        subjectStr = ['Subject' pad(num2str(currSubjId), 2, 'left', '0')];
        sessionStr = ['Session' pad(num2str(currSessId), 2, 'left', '0')];

        currFileEntry.exerciseName = [exerciseNames{ind_exercise}];
        currFileEntry.fileId = databaseSpecTable.(currExerId)(ind_subj);
        currFileEntry.subjectNumber = currSubjId;
        currFileEntry.subjectString = subjectStr;
        currFileEntry.sessionNumber = 1;
        
        
        
        if ~(...
                (currFileEntry.subjectNumber == 13 && strcmpi(currExerId, 'CORK_STD_FAT')))
            
            %         if  ~(currFileEntry.subjectNumber == 03 && strcmpi(currExerId, 'target_55_1_2'))
            % blacklist
            continue
        end
        
%         if currFileEntry.subjectNumber ~= 9
%             continue
%         end
% % %         
%         if ~(strcmpi(currFileEntry.exerciseName, 'SQUA_STD_FAT'))  
% % % % %         if ~(strcmpi(currFileEntry.exerciseName, 'SQUA_STD_FAT')) 
% % %         if ~(strcmpi(currExerId, 'SQUA_STD_CRO') || strcmpi(currExerId, 'SQUA_STD_FAT'))
%             continue
%         end
%         
%         % if the file id is actually empty, skip the entry
%         if isempty(currFileEntry.fileId) || currFileEntry.fileId == 0
%             continue
%         end
        
        %         % rename and copy the mocap databaseSpecTable over to the new destination. if there
        %         % is "trimmed" databaseSpecTable, copy those ones
        %         trcPrefix = fullfile(filepathSource, subjectStr, 'mocap_fp', ['exercise' num2str(currFileEntry.fileId) '_trim*.trc']);
        %         dirTrc = dir(trcPrefix);
        %
        %         % if there is no trimmed trc file, copy the default one instead
        %         if isempty(dirTrc)
        %             trcPrefix = fullfile(filepathSource, subjectStr, 'mocap_fp', ['exercise' num2str(currFileEntry.fileId) '.trc']);
        %             dirTrc = dir(trcPrefix);
        %         end
        trcPrefix = fullfile(filepathSource, subjectStr, 'mocap_fp', ['exercise' num2str(currFileEntry.fileId) '.trc']);
        dirTrc = dir(trcPrefix);
        
        indTrcSuffix = 0;
        for ind_trcFile = 1:length(dirTrc)
            currFile = dirTrc(ind_trcFile).name;
            if length(strsplit(currFile, 'Unnamed')) > 1
                continue
            end
            
            indTrcSuffix = indTrcSuffix + 1;
            [~,currFileName,~] = fileparts(currFile);
            %             targetPath = fullfile(filepathTarget, subjectStr, sessionStr, currFileEntry.exerciseName, 'mocap_fp');
            %             targetPath = fullfile(filepathTarget, 'mocap_fp');
            %             mkdir(targetPath);
            %
            %             % copy over the trc file
            sourceFileTrc = fullfile(filepathSource, subjectStr, 'mocap_fp', [currFileName '.trc']);
            %             targetFileTrc = fullfile(targetPath, [currFileEntry.exerciseName '_' currFileEntry.subjectString '_p' num2str(indTrcSuffix) '.trc']);
            %             copyfile(sourceFileTrc, targetFileTrc);
            %
            %             % now copy over the anc file
            sourceFileAnc = fullfile(filepathSource, subjectStr, 'mocap_fp', [currFileName '.anc']);
            %             targetFileAnc = fullfile(targetPath, [currFileEntry.exerciseName '_' currFileEntry.subjectString '_p' num2str(indTrcSuffix) '.anc']);
            %             copyfile(sourceFileAnc, targetFileAnc);
            
            % run the IK
            currFileEntry.filePathTrc = sourceFileTrc;
            currFileEntry.filePathAnc = sourceFileAnc;
        end
        
        indFileEntry = indFileEntry + 1;
        allFileEntry(indFileEntry) = currFileEntry;
    end
end

%% build model
if runModelParam
    allSubjectString = unique({allFileEntry.subjectString});
    for ind_subj = 1:length(allSubjectString)
        overallModelInd = 0;
        overallModelStruct = [];
        
        currFileEntrySubject = [];
        currFileEntrySubject.subjectString = allSubjectString{ind_subj};
        filepathModelParam = globalConstants_filepaths.fileFullModelParameters(filepathTarget, currFileEntrySubject, externalParam_mocap);
        if exist(filepathModelParam, 'file') && overwriteExisting == 0
            % template already exist, proceed to the next entry
            fprintf('Model template: %s model template at %s already exist, skipping proc\n', currFileEntrySubject.subjectString, filepathModelParam);
            
            if 0
                model = rlModelInstance_iit(0);
                model.loadModelFromModelSpecs(filepathModelXml, filepathModelParam);
                model.plotCurrentPose();
            end
            
            continue
        end
        
        for ind_fileEntry = 1:length(allFileEntry)
            currFileEntry = allFileEntry(ind_fileEntry);
            
            if ~strcmpi(allSubjectString{ind_subj}, currFileEntry.subjectString)
                % if not part of the current subject, ignore
                continue;
            end
            
%             use = blackList_modelParam_iit(filepathBlacklistDecoder, currFileEntry);
%             if ~use
%                 continue;
%             end
            
            fprintf('[%s] runModelMake: Subject %s, Exercise %s\n', datestr(now), currFileEntry.subjectString, currFileEntry.exerciseName);
            
            % generate entry for the model struct
%             try
                externalParam_mocap.ekfRun = 1:20;
                [modelInstance, algorithmParam, featureSet_initPose...
                    dataInstance_trc, dataInstance_ekf, ekfType, ekfTuningParam, initPose, ...
                    kinematicTransform, dynamicTransform, sensorTransform, sensorSecondaryTransform] = initModelData_iit(filepathTarget, currFileEntry, externalParam_mocap);
                
                overallModelInd = overallModelInd + 1;
                if overallModelInd == 1
                    overallModelStruct.kinematicTransform = kinematicTransform;
                    overallModelStruct.dynamicTransform = dynamicTransform;
                    overallModelStruct.sensorTransform = sensorTransform;
                    overallModelStruct.sensorSecondaryTransform = sensorSecondaryTransform;
                else
                    overallModelStruct(overallModelInd).kinematicTransform = kinematicTransform;
                    overallModelStruct(overallModelInd).dynamicTransform = dynamicTransform;
                    overallModelStruct(overallModelInd).sensorTransform = sensorTransform;
                    overallModelStruct(overallModelInd).sensorSecondaryTransform = sensorSecondaryTransform;
                end
%             catch err
%                 err
%                 
%                 filepathModelSpec = globalConstants_filepaths.fileFullKinematicIndividualLog(filepathTarget);
%                 log_modelchar([], filepathModelSpec, currFileEntry, []);
%             end
        end
        
        if ~isempty(overallModelStruct)
            % run the file model production (if appropriate) for the previous
            % subject
            modelStruct = mergeModelData(overallModelStruct);
            save(filepathModelParam, 'modelStruct');
            fprintf('Model template: Saving %s model template at %s\n', currFileEntrySubject.subjectString, filepathModelParam);
            
            % note down what the overall one is
            currFileEntrySubject.exerciseName = 'A_MERGED';
            currFileEntrySubject.fileId = 0;
            filepathModelSpec = globalConstants_filepaths.fileFullKinematicIndividualLog(filepathTarget);
            log_modelchar(modelStruct, filepathModelSpec, currFileEntrySubject, algorithmParam);
            
            [fileID, errmsg] = fopen(filepathModelSpec, 'a');
            fprintf(fileID, '\n');
            fclose(fileID);
        end
    end
end

%% run mocap
if runEkf_mocap
    for ind_fileEntry = 1:length(allFileEntry)
        currFileEntry = allFileEntry(ind_fileEntry);
        fprintf('[%s] runIKEKF: Subject %s, Exercise %s\n', datestr(now), currFileEntry.subjectString, currFileEntry.exerciseName);
        
        try
            externalParam_mocap.ekfRun = [];
            filepathExternalMod = globalConstants_filepaths.fileFullModelParameters(filepathTarget, currFileEntry, externalParam_mocap);
            externalParam_mocap.externalModelDataFilePath = filepathExternalMod;
            main_ik_iit(currFileEntry, filepathTarget, externalParam_mocap);
        catch err
            err
            err.stack
        end
    end
end

if runAnalysis
    currFilepathLogTarget = fullfile(filepathTarget);
    
    % save file names
    filepathMocapIndividualLog = globalConstants_filepaths.fileFullMocapIndividualLog(currFilepathLogTarget);
    filepathMocapGroupLog = globalConstants_filepaths.fileFullMocapGroupLog(currFilepathLogTarget);
    
    % load the model specifics
    id = 0;
    modelInstance = rlModelInstance_iit(id);
    modelInstance.loadModel(filepathModelXml);
    modelInstance.initializeJointAndSensorInds();
    
    mocapTableMaster = [];
    
    for ind_subj = 1:length(databaseSpecTable.SUBJECT)
        currSubjId = databaseSpecTable.SUBJECT(ind_subj);
        
        mocapGroupTableIK = zeros(length(exerciseNames), 1);
        currModelBase = 'floating';
        
        for ind_exercise = 1:length(exerciseNames)
            % load a specific subject and exercise name
            currExerId = exerciseNames{ind_exercise};
            currSessId = 1;
            
            % make sure that it's "subject01" instead of "subject1"
            subjectStr = ['Subject' pad(num2str(currSubjId), 2, 'left', '0')];
            sessionStr = ['Session' pad(num2str(currSessId), 2, 'left', '0')];
            
            % pass into parser
            currFileEntry.exerciseName = [exerciseNames{ind_exercise}];
            currFileEntry.fileId = databaseSpecTable.(currExerId)(ind_subj);
            currFileEntry.subjectNumber = currSubjId;
            currFileEntry.subjectString = subjectStr;
            currFileEntry.sessionNumber = 1;
            
            fprintf('[%s] runAnalysis: Subject %s, Exercise %s\n', datestr(now), subjectStr, currExerId);
            
            % run the fixed base mocap
            currLinkDefinition = 'X00';
            
            % run mocap IK analysis
            try
                markerRmse = calc_error_mocapToTrc_ssa(modelInstance, currLinkDefinition, currFileEntry, ...
                    filepathTarget, filepathMocapIndividualLog, currModelBase);
                mocapGroupTableIK(ind_exercise, 1) = mean(markerRmse);
            catch err
                err
                err.stack
                mocapGroupTableIK(ind_exercise, 1) = -1;
                continue;
            end
        end
        
        localFieldsLinkDef{1} = subjectStr;
        newFieldsY = exerciseNames;
        %             newFieldsY{1, 1} = [currSubjId];
%         newFieldsY{1, 1} = [currExerId];
        T1 = array2table(newFieldsY, 'VariableNames', {'exercise'});
        T2 = array2table(mocapGroupTableIK, 'VariableNames', localFieldsLinkDef);
        if isempty(mocapTableMaster)
            mocapTableMaster = [T1 T2];
        else
            mocapTableMaster = [mocapTableMaster T2];
        end
        
        writetable(mocapTableMaster, filepathMocapGroupLog);
    end
    
    mocapTableMaster
end

