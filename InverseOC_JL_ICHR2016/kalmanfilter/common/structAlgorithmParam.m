function algorithmParam = structAlgorithmParam(structToCheck)
    externalParamParser = inputParser;
    externalParamParser.KeepUnmatched = true;
    
    addOptional(externalParamParser, 'modalityType', 'mocap');
    addOptional(externalParamParser, 'modalityCalibration', 'mocap');
    addOptional(externalParamParser, 'linkDefinition', 'initialpose');
    addOptional(externalParamParser, 'modelBase', 'floating');
    addOptional(externalParamParser, 'subSuffix', '');
    addOptional(externalParamParser, 'datasetId', 'IIT');
    
    addOptional(externalParamParser, 'externalJointAngleData', []);
    addOptional(externalParamParser, 'externalTrcData', []);
    addOptional(externalParamParser, 'externalImuData', []);
    addOptional(externalParamParser, 'externalImuInstance', []);
    addOptional(externalParamParser, 'externalTrcInstance', []);
    addOptional(externalParamParser, 'externalBaseFrameData', []);
    addOptional(externalParamParser, 'initPose', []);
    
    addOptional(externalParamParser, 'externalInitialModelDataFilePath', []);
    addOptional(externalParamParser, 'externalRecoveryModelDataFilePath', []);

    addOptional(externalParamParser, 'linkLengthArrayUse', []);
    addOptional(externalParamParser, 'sensorLengthArrayUse', []);
    addOptional(externalParamParser, 'sensorSecondaryArrayUse', []);
    
    addOptional(externalParamParser, 'ekfTuningParam', 0);
    
    addOptional(externalParamParser, 'overwriteExisting', 0);
    addOptional(externalParamParser, 'visualize_pose', 0);
    addOptional(externalParamParser, 'visualize_main', 0);
    addOptional(externalParamParser, 'endEffectorName', '');
    addOptional(externalParamParser, 'baseFrame', []);
    addOptional(externalParamParser, 'baseMarker', []);
    addOptional(externalParamParser, 'initPoseFrameCount', 100);
    addOptional(externalParamParser, 'missingMarkerValue', 100000);
    addOptional(externalParamParser, 'ekfRun', []);
    addOptional(externalParamParser, 'addUnknownMarkersToEkf', 1);
    addOptional(externalParamParser, 'saveModelLog', 0);
    
    addOptional(externalParamParser, 'imuKinematicChain', []);
    
    if ~exist('structToCheck', 'var')
        structToCheck = struct;  
    end
    
    parse(externalParamParser, structToCheck);
    algorithmParam = externalParamParser.Results;
    
    str = fieldnames(externalParamParser.Unmatched);
    for i = 1:length(str)
        algorithmParam.(str{i}) = externalParamParser.Unmatched.(str{i});
    end
    
    algorithmParam.saveSuffix = globalConstants_filepaths.filePrefixSuffix(algorithmParam);
    
% classdef structAlgorithmParam 
%     properties
%         % algorithm wide params
%         visualize;
%         endEffectorName;
%         ekfRun;
%         ekfForceMatch;
%         
%         baseFrame;
%         baseMarker;
%         
%         % subject specific params
%         subjectId;
%     end
%     
%     methods
% %         function obj = structAlgorithmParam(varargin)
% %             obj.visualize = varargin{1};
% %             obj.endEffectorName = varargin{2};
% %         end
%     end
% end

