% this script loads the content of databaseSpecTablebasespec, runs EKF on the mocap databaseSpecTable, 
% then saves the output into a target folder
clearvars
clc

overwriteExisting = 1;
visualize_pose = 1;
visualize_main = 1;

basepathSource = '.\data\'; 
basepathTarget = '.\results\';

% filepaths to the source databaseSpecTable and target databaseSpecTable
timeStamp = '2018_12_11_singleArm'; filepathModelXml = fullfile('.', 'model', 'ioc_v3_dhAlign.xml');

% timeStamp = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
filepathSource = fullfile(basepathSource,filesep);
filepathTarget = fullfile(basepathTarget,timeStamp,filesep);

% setting the path for common files
addpath(genpath('.\support'));
addpath(genpath('..\common'));
addpath(genpath('..\toolboxes'));
addpath(genpath('..\ekf'));

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

% run inverse kinematics
currFileEntry.exerciseName = 'exercise';
currFileEntry.fileId = [];
currFileEntry.subjectNumber = 1;
currFileEntry.subjectString = '1';
currFileEntry.sessionNumber = 1;
currFileEntry.filePathTrc = '.\data\Subject01_short.trc';
% currFileEntry.filePathTrc = '.\data\Subject01_SingleArm_BaseLine_UBMarkersIncluded.trc';

algorithmParam = structAlgorithmParam(externalParam_mocap);

filepathModelInitPose = filepathModelXml;
filepathModelRecoveryTRC = filepathModelXml;

baseFrameOrig = 'frame_6dof_root';
baseFrameCurr = 'frame_6dof_root'; % frame_6dof_root frame_6dof2_root
baseFrameMarker = ''; % HIP_JC_R ANKLE_JC_R

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
modelInstance = rlModelInstance_ioc(currFileEntry.subjectNumber);
dataInstance_trc = rlDataInstance_trc(modelInstance);
dataInstance_trc.loadData(currFileEntry.filePathTrc, algorithmParam);
dataInstance_trc.data.SHOULDER_R_JC = ones(size(dataInstance_trc.data.SHOULDER_R_JC))*1e-6; % generally speaking, the data does not like zeros
dataInstance_trc.dataProcessing(algorithmParam);

% link the data instance into the model and create model to init XML model
% linkage length and sensor attachments
[kinematicTransform, dynamicTransform, sensorTransform] = ...
    modelInstance.makeModel(filepathModelInitPose, dataInstance_trc, algorithmParam);

mdl = modelInstance.model;
mdl.forwardPosition();

vis = rlVisualizer('vis',640,480);
mdl.forwardPosition();
vis.addModel(mdl);
applyMarkersToVisualization(vis, mdl, [], [], [], []);
vis.update();
  
mdl.position = [0 0 0 0 0];
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mdl.inverseDynamics();
mdl.torque'

mdl.position = [0 pi/2 0 0 0];
mdl.forwardPosition();
mdl.forwardVelocity();
mdl.forwardAcceleration();
mdl.inverseDynamics();
mdl.torque'