function [modelInstance, algorithmParam, featureSet_initPose...
    dataInstance_trc, dataInstance_ekf, ekfType, ekfTuningParam, initPose, ...
    kinematicTransform, dynamicTransform, sensorTransform, sensorSecondaryTransform] = initModelData_expressiveioc(filepathTarget, currFileEntry, pExternalParam, filepathModelXml)

    algorithmParam = structAlgorithmParam(pExternalParam);

    ekfTypeInitPose = 'GP_EKF_Q_DQ';
    switch algorithmParam.modalityType
        case 'mocap'
            ekfType = 'GP_EKF_Q_DQ';
    end
    
    switch algorithmParam.modelBase
        case 'floating'
            filepathModelInitPose = filepathModelXml;
            filepathModelRecoveryTRC = filepathModelXml;

            baseFrameOrig = 'frame_6dof_root';
            baseFrameCurr = 'frame_6dof_root'; % frame_6dof_root frame_6dof2_root
            baseFrameMarker = ''; % HIP_JC_R ANKLE_JC_R
    end

    % algorithmParam = structAlgorithmParam;
    algorithmParam.endEffectorName = '';
    algorithmParam.baseFrameOrig = baseFrameOrig;
    algorithmParam.baseFrame = baseFrameCurr;
    algorithmParam.baseMarker = baseFrameMarker; 
    algorithmParam.initPoseFrameCount = 100;
    algorithmParam.addUnknownMarkersToEkf = 1;
    % algorithmParam.ekfRun = 1:200;%500:1500;
    algorithmParam.subjectId = currFileEntry.subjectNumber;

    % initialize the model instance and load marker data
    modelInstance = rlModelInstance_expressiveioc(currFileEntry.subjectNumber);
    % modelInstance = rlModelInstance_expressiveioc_rightArm(currFileEntry.subjectNumber);
    dataInstance_trc = rlDataInstance_trc(modelInstance);
    dataInstance_trc.loadData(currFileEntry.filePathTrc, algorithmParam);
%     dataInstance_trc.data.SHOULDER_R_JC = ones(size(dataInstance_trc.data.SHOULDER_R_JC))*1e-6; % generally speaking, the data does not like zeros
    dataInstance_trc.dataProcessing(algorithmParam);

    % link the data instance into the model and create model to init XML model
    % linkage length and sensor attachments
    [kinematicTransform, dynamicTransform, sensorTransform] = ...
        modelInstance.makeModel(filepathModelInitPose, dataInstance_trc, algorithmParam);

%     modelInstance.initializeJointAndSensorInds();
    modelInstance.plotCurrentPose();
    
    modelInstance.model.base = algorithmParam.baseFrame;
    dataInstance_trc.updateBaseMarker(algorithmParam.baseMarker);
    
    %% inititalize the EKF model
    ekfTuningParam.processNoiseCoefficient = 1e6;
    ekfTuningParam.processNoiseCoefficientPosition = 0;
    ekfTuningParam.processNoiseCoefficientVelocity = 0;
    ekfTuningParam.processNoiseCoefficientAcceleration = 0;
    ekfTuningParam.observationNoiseCoefficientMarkerPosition = 0.01;
    ekfTuningParam.observationNoiseCoefficientMarkerVelocity = 0;
    ekfTuningParam.observationNoiseCoefficientAccelerometer = 0;
    ekfTuningParam.observationNoiseCoefficientGyroscope = 0;
    ekfTuningParam.observationNoiseCoefficientYaw = 0;
    ekfTuningParam.markerSwapTolerance = 100;
    
    % initialize the EKF state by setting an initial pose for the 3 DOF floating
    % base by feeding the first frame into EKF 100 times
    indFrameToRepeat = 1;
    [timeArray, dataArray, indArray, labelArray, matchArray] = ...
        dataInstance_trc.repeatMeasurementForEkf(modelInstance.model, indFrameToRepeat, algorithmParam.initPoseFrameCount); % repeat the first frame 100 times
    % matchArray = dataInstance_trc.outputMatchMatrixCheckDataInstance(modelInstance.model, modelInstance);
    
    featureSet_initPose = rlFeatureSet_ioc();
    featureSet_initPose.setModelParam(modelInstance.model, dataInstance_trc.dt, algorithmParam, modelInstance);
    featureSet_initPose.setEKFParams(modelInstance, ekfTypeInitPose, ekfTuningParam);
    featureSet_initPose.algorithmParam.ekfForceMatch = 1;
    featureSet_initPose.algorithmParam.visualize = algorithmParam.visualize_pose;
    featureSet_initPose.baseFrame = 'frame_6dof_root';
    
    initPos = zeros(size(modelInstance.model.position));
%     initPosHip = dataInstance_trc.data.SHOULDER_R_JC(1, :);
%     if sum(find(isnan(initPosHip) == 1)) == 0
%         initPos(1:3) = initPosHip;
%     end
%     modelInstance.model.position = initPos;
    
    featureSet_initPose.calcStatesFromTrc_BaseFrame(dataArray, initPos, dataInstance_trc, modelInstance, [], indArray, labelArray, matchArray);
    initPose = featureSet_initPose.model.position;
    jointAngles = repmat(initPose', length(modelInstance.lengthUseInd), 1);
    
    % set up model parameters by copying from init frame model
    modelInstance.model.base = baseFrameOrig;
    sensorSecondaryTransform = modelInstance.makeTRCModel(filepathModelRecoveryTRC, dataInstance_trc, initPose, jointAngles, algorithmParam, []);
    modelInstance.model_ekf.base = baseFrameCurr;
    
    % TRC data already loaded in dataInstance
    dataInstance_ekf = dataInstance_trc;
end