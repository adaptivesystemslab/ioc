function currModelStruct = modelSpecFromEkfMatches(basepathTarget, currFileEntry, featureSet, pExternalParam)
    % given the featureset match matrix and dataset, produce an updated
    % kinematic transform set 
    
    algorithmParam = structAlgorithmParam(pExternalParam);
    algorithmParam.saveModelLog = 1;
             
    filepathModelParamStruct = globalConstants_filepaths.fileFullModelParameterStruct(basepathTarget, currFileEntry, algorithmParam);  
    
    if exist(filepathModelParamStruct, 'file') && algorithmParam.overwriteExisting == 0
        load(filepathModelParamStruct);
        return;
    end

    switch algorithmParam.datasetId
        case 'IIT'
            modelInstance = rlModelInstance_iit(currFileEntry.subjectNumber);
            filepathModelRecoveryIMU = fullfile('.', 'model', 'iit_v10.xml');

        case 'SSA'
            modelInstance = rlModelInstance_ssa(currFileEntry.subjectNumber);
            filepathModelRecoveryIMU = fullfile('.', 'model', 'ssa_v4_floating.xml');
            
        case 'GAIT'
            modelInstance = rlModelInstance_gait(currFileEntry.subjectNumber);
            filepathModelRecoveryIMU = fullfile('.', 'model', 'avatar_v04_p_bothLegs_zyx_hip.xml');
    end
    
    
%     filepathModelInitPose = fullfile('.', 'model', 'ssa_v4_floating.xml');
%     filepathModelRecoveryTRC = fullfile('.', 'model', 'ssa_v4_floating.xml');
    
    
%     dataInstance_trc = rlDataInstance_trc(modelInstance);
%     
%     dataInstance_trc.loadData(currFileEntry.filePathTrc, algorithmParam);
%     dataInstance_trc.dataProcessing(algorithmParam);
    
    % determining which matrix matches
    [ekf_markerMask, ekf_eventhandler, ekf_markerTally] = ekfSensorMatching(featureSet);
    correctArray = ekf_markerMask.correct;
    
    masterLengthIndsUse = modelInstance.lengthUseInd;
    if masterLengthIndsUse == 0
        masterLengthIndsUse = 1:size(correctArray, 1);
    end
    
    % generate in linkLengthArray sequence the indices that can be used
    % iterate though the match matrix with the indices that it's suppose to be    
    % link length composition markers
    linkLengthArrayUse = [];
    for i = 1:length(modelInstance.linkLengthArray)
        sourceInds = findUsableIndsGroup(correctArray, modelInstance.linkLengthArray{i}{3}, modelInstance, featureSet);
        targetInds = findUsableIndsGroup(correctArray, modelInstance.linkLengthArray{i}{4}, modelInstance, featureSet);
        
        useInds1 = intersect(sourceInds, targetInds);
        useInds2 = intersect(useInds1, masterLengthIndsUse);
        linkLengthArrayUse{i} = useInds2;
    end
    
    % sensor composition markers
    sensorLengthArrayUse = [];
    for i = 1:length(modelInstance.flatSensorAttachmentArray)
        sourceInds = findUsableIndsGroup(correctArray, modelInstance.flatSensorAttachmentArray{i}{2}, modelInstance, featureSet);
        targetInds = findUsableIndsSingle(correctArray, modelInstance.flatSensorAttachmentArray{i}{3}, modelInstance, featureSet);
        
        useInds1 = intersect(sourceInds, targetInds);
        useInds2 = intersect(useInds1, masterLengthIndsUse);
        sensorLengthArrayUse{i} = useInds2;
    end
    
    for i = 1:length(modelInstance.flatSensorSecondaryAttachmentArray)
        flatSensorSecondaryLastCell(i) = modelInstance.flatSensorSecondaryAttachmentArray{i}(end);
    end
    
    % secondary sensors composition markers
    sensorSecondaryArrayUse = [];
    for i = 1:length(modelInstance.sensorSecondaryAttachmentArray)
        sourceInds = findUsableIndsGroup(correctArray, modelInstance.sensorSecondaryAttachmentArray{i}{2}, modelInstance, featureSet);
        
        switch numel(modelInstance.sensorSecondaryAttachmentArray{i})
            case {3}
                targetInds1 = findUsableIndsSingle(correctArray, modelInstance.sensorSecondaryAttachmentArray{i}{3}, modelInstance, featureSet);
                indsFinal = intersect(sourceInds, targetInds1);
                
            case {5, 6}
                targetInds1 = findUsableIndsSingle(correctArray, modelInstance.sensorSecondaryAttachmentArray{i}{3}, modelInstance, featureSet);
                targetInds2 = findUsableIndsSingle(correctArray, modelInstance.sensorSecondaryAttachmentArray{i}{4}, modelInstance, featureSet);
                targetInds3 = findUsableIndsSingle(correctArray, modelInstance.sensorSecondaryAttachmentArray{i}{5}, modelInstance, featureSet);
                inds1 = intersect(sourceInds, targetInds1);
                inds2 = intersect(inds1, targetInds2);
                inds3 = intersect(inds2, targetInds3);
                indsFinal = intersect(inds3, masterLengthIndsUse);
        end
        
%         for k = 3:length(modelInstance.sensorSecondaryAttachmentArray{i})
%             indFlatSensor = find(ismember(flatSensorSecondaryLastCell, modelInstance.sensorSecondaryAttachmentArray{i}{k}));
%             sensorSecondaryArrayUse{indFlatSensor} = indsFinal;
%         end
        
        sensorSecondaryArrayUse{i} = indsFinal;
    end
    
    % with the new induse arrays built, pass it back into initModel for
    % processing
    algorithmParam.linkLengthArrayUse = linkLengthArrayUse;
    algorithmParam.sensorLengthArrayUse = sensorLengthArrayUse;
    algorithmParam.sensorSecondaryArrayUse = sensorSecondaryArrayUse;
    
    switch algorithmParam.datasetId
        case 'IIT'
            [modelInstance, algorithmParam, dataInstance_trc, dataInstance_ekf, qInit, dqInit, ddqInit, baseFrameCurr, fullPose, ...
                kinematicTransform, dynamicTransform, sensorTransform, sensorSecondaryTransformTrc] = ...
                initModelData_iit(basepathTarget, currFileEntry, algorithmParam);
            
            sensorSecondaryTransformImu = [];
            
        case 'SSA'
            [modelInstance, algorithmParam, dataInstance_trc, dataInstance_ekf, qInit, dqInit, ddqInit, baseFrameCurr, fullPose, ...
                kinematicTransform, dynamicTransform, sensorTransform, sensorSecondaryTransformTrc] = ...
                initModelData_ssa(basepathTarget, currFileEntry, algorithmParam);
            
            R_sensor = eye(3);
            sensorSecondaryTransformImu = modelInstance.makeIMUModel(filepathModelRecoveryIMU, ...
                dataInstance_trc, R_sensor, qInit, fullPose, algorithmParam, sensorSecondaryArrayUse);
            
        case 'GAIT'
            [modelInstance, algorithmParam, dataInstance_trc, dataInstance_ekf, qInit, dqInit, ddqInit, baseFrameCurr, fullPose, ...
                kinematicTransform, dynamicTransform, sensorTransform, sensorSecondaryTransformTrc] = ...
                initModelData_gait(basepathTarget, currFileEntry, algorithmParam);
            
            R_sensor = eye(3);
            sensorSecondaryTransformImu = modelInstance.makeIMUModel(filepathModelRecoveryIMU, ...
                dataInstance_trc, R_sensor, qInit, fullPose, algorithmParam, sensorSecondaryArrayUse);            
    end
   
    clear featureSet;
    clear dataInstance_trc;
    clear modelInstance;
    
    currModelStruct.exerciseName = currFileEntry.exerciseName;
    currModelStruct.kinematicTransform = kinematicTransform;
    currModelStruct.dynamicTransform = dynamicTransform;
    currModelStruct.sensorTransform = sensorTransform;
    currModelStruct.sensorSecondaryTransformTrc = sensorSecondaryTransformTrc;
    currModelStruct.sensorSecondaryTransformImu = sensorSecondaryTransformImu;
    
    save(filepathModelParamStruct, 'currModelStruct');

%     % pass it into makemodel and sensortransform
%     % update kin/dyn/sen trans
%     [kinematicTransform, dynamicTransform, sensorTransform] = ...
%             modelInstance.makeModel(filepathModelInitPose, dataInstance_trc, algorithmParam, linkLengthArrayUse, sensorLengthArrayUse);
%     
%     % apply the first posture to the model
%     initPose = featureSet.q(1, :);
%     modelInstance.model.position = featureSet.q(1, :);
%     modelInstance.model.forwardPosition();
%     
%      sensorSecondaryTransform = modelInstance.makeTRCModel(filepathModelRecoveryTRC, ...
%                     dataInstance_trc, initPose, jointAngles, algorithmParam);
%     % calculate 
% %       switch algorithmParam.modalityType
% %             case 'imu'
% %                 sensorSecondaryTransform = modelInstance.makeIMUModel(filepathModelRecoveryIMU, ...
% %                     dataInstance_trc, R_sensor, initPose, jointAngles, algorithmParam);
% %                 modelInstance.model_ekf.base = baseFrameCurr;
% % 
% %             case 'mocap'
% %                 sensorSecondaryTransform = modelInstance.makeTRCModel(filepathModelRecoveryTRC, ...
% %                     dataInstance_trc, initPose, jointAngles, algorithmParam);
% %         end
% 
%     
%     % (also make an IMU version so we don't need to run it twice)
%     
    
end