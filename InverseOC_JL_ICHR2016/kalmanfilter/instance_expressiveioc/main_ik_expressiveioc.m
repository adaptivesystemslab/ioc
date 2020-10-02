function main_ik_expressiveioc(currFileEntry, filepathInstance, pExternalParam, filepathModelXml)
%% IK EKF code for IIT dataset
algorithmParam = structAlgorithmParam(pExternalParam);

% if file exist already, don't run it
savePath = globalConstants_filepaths.fileFullEkfIk(filepathInstance, currFileEntry, algorithmParam);
if exist(savePath, 'file') && ~algorithmParam.overwriteExisting
    return
end

[modelInstance, algorithmParam, featureSet_initPose...
    dataInstance_trc, dataInstance_ekf, ekfType, ekfTuningParam, initPose] = initModelData_expressiveioc(filepathInstance, currFileEntry, algorithmParam, filepathModelXml);

%% run EKF
[timeArray, dataArray, indArray, labelArray, matchArray] = dataInstance_ekf.outputMeasurementForEkf(modelInstance.model_ekf);

featureSet_recovery = rlFeatureSet_ioc();
featureSet_recovery.visSkipRate = 1;
featureSet_recovery.setModelParam(modelInstance.model_ekf, dataInstance_ekf.dt, algorithmParam, modelInstance);
featureSet_recovery.setEKFParams(modelInstance, ekfType, ekfTuningParam);
featureSet_recovery.algorithmParam.visualize = algorithmParam.visualize_main;
featureSet_recovery.algorithmParam.ekfForceMatch = 0; % -1: force match only first frame

frameNames = {'body_base', ...
    'frame_l5_4', 'frame_t10_4', 'frame_c7_6', 'frame_c1head_end', ...
    'frame_c7rshoulder_2', 'frame_rshoulder_6', 'frame_relbow_4', 'frame_rwrist_4', 'frame_rhand_end', ...
    'frame_c7lshoulder_2', 'frame_lshoulder_6', 'frame_lelbow_4', 'frame_lwrist_4', 'frame_lhand_end'};
% frameNames = {'body_base', 'frame_rshoulder_6', 'frame_relbow_4', 'frame_rwrist_4', 'frame_rhand_end'};
featureSet_recovery.setFrameSaving(frameNames, size(dataArray, 1));

featureSet_recovery.calcStatesFromTrc_BaseFrame(dataArray, initPose, dataInstance_trc, modelInstance, algorithmParam, indArray, labelArray, matchArray);

%% Post-processing after EKF
numJoints = length(modelInstance.model_ekf.joints); % removing 6DOF free joint at the start of the model
inds_fb = 0;
inds_fb_pos = 1:inds_fb;
inds_model_pos = (inds_fb+1):numJoints;
inds_fb_vel = (numJoints+1):(numJoints+inds_fb);
inds_model_vel = (numJoints+inds_fb+1):(2*numJoints);

featureSet_recovery.calcQfromEkfStates(inds_fb_pos, inds_model_pos, inds_fb_vel, inds_model_vel);
featureSet_recovery.calcFeaturesFromQ(timeArray);

% rearrange measurement_input and measurement_output so it matches the
% intended array from rlmodelInstance_iit.sensorNameFull before saving
% it to file
measurement_label_new = modelInstance.sensorNameFull;
[measurement_input_new, input_matching] = modelInstance.rearrangeMeasurement(dataArray, labelArray, measurement_label_new);
[measurement_output_new, output_matching] = modelInstance.rearrangeMeasurement(featureSet_recovery.measurement_output, featureSet_recovery.measurement_labels, ...
    measurement_label_new, measurement_input_new);

featureSet_recovery.measurement_labels = measurement_label_new;
featureSet_recovery.measurement_input = measurement_input_new;
featureSet_recovery.measurement_output = measurement_output_new;
featureSet_recovery.measurement_input_match = input_matching;
featureSet_recovery.measurement_output_match = output_matching;
featureSet_recovery.measurement_ekf_match = matchArray(:, 2);

[ekf_markerMask, ekf_eventhandler, ekf_markerTally] = ekfSensorMatching(featureSet_recovery);
featureSet_recovery.measurement_mask = ekf_markerMask.swappedmissing;

if ~isempty(filepathInstance)
    fprintf('File saved to %s\n', savePath);
    savePath = globalConstants_filepaths.fileFullEkfIk(filepathInstance, currFileEntry, algorithmParam);
    featureSet_recovery.segments = dataInstance_trc.segments;
    saveVar = featureSet_recovery.saveData(savePath, algorithmParam.baseFrame, ekfTuningParam, modelInstance);
end

% clear objects
clear featureSet_initPose featureSet_recovery;
clear dataInstance_clean dataInstance_ekf dataInstance_gt;
clear modelInstance;