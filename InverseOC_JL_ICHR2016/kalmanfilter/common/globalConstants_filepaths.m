classdef globalConstants_filepaths 
    % create a central method to control the filepath output
    properties(Constant = true)
        % directory prefixes
        folderSuffixCalibration = 'matCalibration';
        folderSuffixTimeAlignment = 'figTimeAlignment';
    end
    
    methods(Static = true)              
        %
        function [folderPath] = folderCalibrationBase(basePath)
            folderPath = fullfile(basePath, globalConstants_filepaths.folderSuffixCalibration);
            checkMkdir(folderPath);
        end
        
        function [folderPath] = folderLogBase(basePath)
            folderPath = fullfile(basePath, 'log');
            checkMkdir(folderPath);
        end
        
        %%
        function [filePath, folderPath] = fileTimeAlignedLog(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, globalConstants_filepaths.folderSuffixTimeAlignment);
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_timeAlignmentLog.mat'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullEkfIk(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'matEkfIk', globalConstants_filepaths.filePrefixSuffix(algorithmParam));
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_ekfIk.mat'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullEkfFk(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'matEkfFk');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_ekfFk.mat'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullModelParametersTrc(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'modelParam');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixSubject(currFileEntry, algorithmParam);
            fileName = [filePrefix '_mdlParamTrc.mat'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullModelParameterStruct(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'modelStruct');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_mdlStruct.mat'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullModelParametersImu(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'modelParam');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixSubject(currFileEntry, algorithmParam);
            fileName = [filePrefix '_mdlParamImu.mat'];
            filePath = fullfile(folderPath, fileName);
        end        
        
        function [filePath, folderPath] = fileFullModelParametersExerciseSpecificTrc(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'modelParam');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixSubjectSpecific(currFileEntry, algorithmParam);
            fileName = [filePrefix '_mdlParamTrc.mat'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullModelParametersExerciseSpecificImu(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'modelParam');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixSubjectSpecific(currFileEntry, algorithmParam);
            fileName = [filePrefix '_mdlParamImu.mat'];
            filePath = fullfile(folderPath, fileName);
        end
        
        %%
        function [filePath, folderPath] = fileFullKinematicIndividualLog(basePath)
            folderPath = fullfile(basePath, 'log');
            checkMkdir(folderPath);
            
            fileName = ['kinematic_mocap_individual.csv'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullMocapIndividualLog(basePath, suffix)
            folderPath = fullfile(basePath, 'log');
            checkMkdir(folderPath);
            
            if ~exist('suffix', 'var')
                suffix = '';
            end
            
            fileName = ['metrics_mocap_individual' suffix '.csv'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullMocapMocapIndividualLog(basePath, suffix)
            folderPath = fullfile(basePath, 'log');
            checkMkdir(folderPath);
            
            if ~exist('suffix', 'var')
                suffix = '';
            end
            
            fileName = ['metrics_mocap_mocap_individual' suffix '.csv'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullMocapGroupLog(basePath, suffix)
            folderPath = fullfile(basePath, 'log');
            checkMkdir(folderPath);
            
            if ~exist('suffix', 'var')
                suffix = '';
            end
            
            fileName = ['metrics_mocap_group' suffix '.csv'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullImuIndividualLog(basePath)
            folderPath = fullfile(basePath, 'log');
            checkMkdir(folderPath);

            fileName = ['metrics_imu_individual.csv'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullImuGroupLog(basePath)
            folderPath = fullfile(basePath, 'log');
            checkMkdir(folderPath);
            
            fileName = ['metrics_imu_group.csv'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullImuMocapFkIndividualLog(basePath)
            folderPath = fullfile(basePath, 'log');
            checkMkdir(folderPath);
            
            fileName = ['metrics_imufkmocap_individual.csv'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullImuMocapFkGroupLog(basePath)
            folderPath = fullfile(basePath, 'log');
            checkMkdir(folderPath);
            
            fileName = ['metrics_imufkmocap_group.csv'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullFileExistLog(basePath)
            folderPath = fullfile(basePath, 'log');
            checkMkdir(folderPath);
            
            fileName = ['metrics_file_exists.csv'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = fileFullFileRuntimeLog(basePath)
            folderPath = fullfile(basePath, 'log');
            checkMkdir(folderPath);
            
            fileName = ['runtimeLog.csv'];
            filePath = fullfile(folderPath, fileName);
        end

        %%
        function [filePath, folderPath] = filePrefixErrorCheckPlot(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'figErrorCheck');
            checkMkdir(folderPath);
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_errorCheck'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = filePrefixErrorCheckPlot2(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'figErrorCheck2');
            checkMkdir(folderPath);
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_errorCheck'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = filePrefixErrorCheckPlot3(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'figErrorCheck3');
            checkMkdir(folderPath);
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_errorCheck'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = filePrefixMatchMatrix(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'figEkfMatchMatrix');
            checkMkdir(folderPath);
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_matchMatrix'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = filePrefixJointsCombined(basePath, currFileEntry, algorithmParam)
            folderPrefix = globalConstants_filepaths.filePrefixTypeSpecific(currFileEntry, algorithmParam);
            
            folderPath = fullfile(basePath, 'figEkfJoints', folderPrefix, 'combined_q');
            checkMkdir(folderPath);
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_jointsCombined'];
            filePath = fullfile(folderPath, fileName);
        end   
        
        function [filePath, folderPath] = filePrefixVelocitiesCombined(basePath, currFileEntry, algorithmParam)
            folderPrefix = globalConstants_filepaths.filePrefixTypeSpecific(currFileEntry, algorithmParam);
            
            folderPath = fullfile(basePath, 'figEkfJoints', folderPrefix, 'combined_dq');
            checkMkdir(folderPath);
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_velocitiesCombined'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = filePrefixJointsIndividual(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'figEkfJoints', 'individual');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_jointsIndividual'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = filePrefixMeasurementsIndividual(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'figEkfMeas');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_mes'];
            filePath = fullfile(folderPath, fileName);
        end        
        
        
        function [filePath, folderPath] = filePrefixMocapImuIkRmse(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'figIkRmseComb');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_ikRmse'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = filePrefixMocapImuIkRmseIndividual(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'figIkRmseInd');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_ikRmse'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = filePrefixMocapImuTimeAlignment(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'figTimeAlignment');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_timeAlign'];
            filePath = fullfile(folderPath, fileName);
        end     
        
        function [filePath, folderPath] = filePrefixFkIndividual(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'figFkMeas');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_mes'];
            filePath = fullfile(folderPath, fileName);
        end

        function [filePath, folderPath] = filePrefixRotDistIndividual(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'figIkRot');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_rotPos'];
            filePath = fullfile(folderPath, fileName);
        end
        
        function [filePath, folderPath] = filePrefixRotVelIndividual(basePath, currFileEntry, algorithmParam)
            folderPath = fullfile(basePath, 'figIkRot');
            checkMkdir(folderPath);
            
            filePrefix = globalConstants_filepaths.filePrefixCombined(currFileEntry, algorithmParam);
            fileName = [filePrefix '_rotVel'];
            filePath = fullfile(folderPath, fileName);
        end
        
        %%
        function fileSuffix = filePrefixSuffix(algorithmParam)
            if ~isempty(algorithmParam.subSuffix)
                fileSuffix = [algorithmParam.modalityType ...
                    '_' algorithmParam.modalityCalibration ...
                    '_' algorithmParam.linkDefinition ...
                    '_' algorithmParam.modelBase ...
                    '_' algorithmParam.subSuffix];
            else
                fileSuffix = [algorithmParam.modalityType ...
                    '_' algorithmParam.modalityCalibration ...
                    '_' algorithmParam.linkDefinition ...
                    '_' algorithmParam.modelBase];
            end
            
        end
        
        function fileSuffix = filePrefixSubjectSpecific(currFileEntry)
            fileSuffix = [currFileEntry.exerciseName '_' currFileEntry.subjectString];
        end
        
        function fileSuffix = filePrefixTypeSpecific(currFileEntry, algorithmParam)
            if isfield(algorithmParam, 'ekfTuningParam')
                if isfield(algorithmParam.ekfTuningParam, 'id')
                    ekfId = algorithmParam.ekfTuningParam.id;
                else
                    ekfId = 1;
                end
            else
                ekfId = 1;
            end
            
            fileSuffix = [globalConstants_filepaths.filePrefixSuffix(algorithmParam) '_' ...
                'ekfId' num2str(ekfId)];
        end
        
        function fileSuffix = filePrefixCombined(currFileEntry, algorithmParam)
            if isfield(algorithmParam, 'ekfTuningParam')
                if isfield(algorithmParam.ekfTuningParam, 'id')
                    ekfId = algorithmParam.ekfTuningParam.id;
                else
                    ekfId = 1;
                end
            else
                ekfId = 1;
            end
            
%             fileSuffix = [globalConstants_filepaths.filePrefixSubjectSpecific(currFileEntry) '_' ...
%                 globalConstants_filepaths.filePrefixSuffix(algorithmParam) '_' ...
%                 'ekfId' num2str(ekfId)];

            fileSuffix = [globalConstants_filepaths.filePrefixSubjectSpecific(currFileEntry) '_' ...
               currFileEntry.sessionString];
        end
        
        function fileSuffix = filePrefixSubject(currFileEntry, algorithmParam)
            fileSuffix = [currFileEntry.subjectString '_' ...
                globalConstants_filepaths.filePrefixSuffix(algorithmParam)];
        end
    end
end