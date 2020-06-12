% this script loads the content of databaseSpecTablebasespec, runs EKF on the mocap databaseSpecTable, 
% then saves the output into a target folder
clearvars
clc

runModelParam = 0;
runEkf_mocap = 1;
runAnalysis = 1;

overwriteExisting = 0;
visualize_pose = 0;
visualize_main = 0;

% basepathSource = 'D:\aslab\data\PamelasData'; 
% basepathTarget = 'D:\aslab\data_IK\PamelasData';
basepathSource = '../../../data/Vicon'; 
basepathTarget = '../../../data/IK_BothArmsCondition';

% filepaths to the source databaseSpecTable and target databaseSpecTable
% timeStamp = '2019_04_11_fullbody3'; filepathModelXml = fullfile('.', 'model', 'ioc_v4_upperbody.xml'); folderExt = '\Full_body';
% timeStamp = '2019_05_27_fullbody'; filepathModelXml = fullfile('.', 'model', 'ioc_v4_rightarm.xml'); folderExt = '/Subject01/BothArms';
timeStamp = '2019_05_27_fullbody'; filepathModelXml = fullfile('', 'model', 'ioc_v4_upperbody.xml'); folderExt = '/Subject01/BothArms';

% timeStamp = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
filepathSource = fullfile(basepathSource,filesep);
filepathTarget = fullfile(basepathTarget,timeStamp,filesep);

% setting the path for common files
addpath(genpath('../toolboxes'));
addpath(genpath('./support'));
addpath(genpath('../common'));
addpath(genpath('../../General_FKEKF/DynamicsModelMatlab/MatlabWrapper'));

externalParam_mocap = [];
externalParam_mocap.modalityType = 'mocap'; % mocap imu
externalParam_mocap.modalityCalibration = 'mocap'; % setting to set calibration matrices
externalParam_mocap.linkDefinition = 'X00'; % setting to set kinematic for trc model (initialpose X00 anth)
externalParam_mocap.modelBase = 'floating'; % floating hip
externalParam_mocap.subSuffix = '';
externalParam_mocap.externalTrcData = [];
externalParam_mocap.externalImuData = [];
externalParam_mocap.externalBaseFrameData = [];
externalParam_mocap.overwriteExisting = overwriteExisting;
externalParam_mocap.visualize_pose = visualize_pose;
externalParam_mocap.visualize_main = visualize_main;

% run inverse kinematics
indFileEntry = 0;
dirSource = dir(filepathSource);
for i = 1:length(dirSource)
    if ~dirSource(i).isdir
        continue
    end
    
    dirSubDir = dir([filepathSource dirSource(i).name folderExt]);  
    for j = 1:length(dirSubDir)
        [filepath,name,ext] = fileparts(dirSubDir(j).name);
        if strcmpi(ext, '.trc')
            currFileEntry.exerciseName = name;
            currFileEntry.fileId = [];
            currFileEntry.subjectNumber = str2num(name(strfind(name,'_')-2:strfind(name,'_')-1)); %str2num(dirSource(i).name(end-1:end));
            currFileEntry.subjectString = name(strfind(name,'Subject'):strfind(name,'_')-1);%dirSource(i).name(end-1:end);
            currFileEntry.sessionNumber = 1;
            % currFileEntry.filePathTrc = '.\data\Subject01SingleArm_BaseLineModified_UBMarkers.trc';
            currFileEntry.filePathTrc = [dirSubDir(j).folder '/' dirSubDir(j).name];
            
            indFileEntry = indFileEntry + 1;
            allFileEntry(indFileEntry) = currFileEntry;
        end
    end
end

for ind_fileEntry = 1:length(allFileEntry)
    currFileEntry = allFileEntry(ind_fileEntry);
    fprintf('[%s] runIKEKF: Subject %s, Exercise %s\n', datestr(now), currFileEntry.subjectString, currFileEntry.exerciseName);
    
    externalParam_mocap.ekfRun = [];
    try
        main_ik_expressiveioc(currFileEntry, filepathTarget, externalParam_mocap, filepathModelXml);
    catch err
        err
    end
end

% modelInstance = rlModelInstance_expressiveioc(0);
modelInstance = rlModelInstance_expressiveioc_rightArm(0);
modelInstance.loadModel(filepathModelXml);
for ind_fileEntry = 1:length(allFileEntry)
    currFileEntry = allFileEntry(ind_fileEntry);
    
    [featureSet_mocap, pathToSave_featureSetEkf_ik, pathToSave_featureSetEkf_fk] = ...
        loadFSOS(currFileEntry, filepathTarget, externalParam_mocap);
    [ekf_markerMask, ekf_eventhandler, ekf_markerTally] = ekfSensorMatching(featureSet_mocap);
     a
    % plot the outputs
    plot_batch_featureSet(currFileEntry, filepathTarget, featureSet_mocap, ...
        ekf_markerMask, ekf_eventhandler, ekf_markerTally, modelInstance, externalParam_mocap);
end

