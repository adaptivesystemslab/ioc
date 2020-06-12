function [initPose, featureSet_initPose] = ekfInitPose(modelInstance, mdl, dataInstance, algorithmParam)
    indFrameToRepeat = 1;
    ekfTypeInitPose = 'GP_EKF_Q_DQ'; % init pose always from mocap
    
    [timeArray, dataArray, indArray, labelArray, matchArray] = ...
        dataInstance.repeatMeasurementForEkf(mdl, indFrameToRepeat, algorithmParam.initPoseFrameCount); % repeat the first frame 100 times
    % [matchMatrix] = dataInstance_trc.outputMatchMatrix(modelInstance.model);

    featureSet_initPose = rlFeatureSet_ioc();
    featureSet_initPose.setModelParam(mdl, dataInstance.dt, algorithmParam, modelInstance);
    featureSet_initPose.setEKFParams(modelInstance, ekfTypeInitPose, algorithmParam.ekfTuningParam);
    featureSet_initPose.algorithmParam.ekfForceMatch = -2;
    featureSet_initPose.algorithmParam.visualize = algorithmParam.visualize_pose;
%     featureSet_initPose.baseFrame = 'body_base';

    featureSet_initPose.calcStatesFromTrc_BaseFrame(dataArray, [], [], [], dataInstance, modelInstance, indArray, labelArray, matchArray);
    initPose = featureSet_initPose.model.position';
end