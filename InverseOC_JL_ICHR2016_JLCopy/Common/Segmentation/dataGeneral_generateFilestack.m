function [fileStack, exerciseStruct] = dataGeneral_generateFilestack(pathToRawData, dataPackageUncheck)
% given the datapackage, this function pulls out all the file information
% that matches the outlined in the datapackage. 

fileStack = {};

if ~iscell(dataPackageUncheck)
    dataPackage{1} = dataPackageUncheck;
else
    dataPackage = dataPackageUncheck;
end

for ind_package = 1:length(dataPackage)
    specStruct = dataPackage{ind_package};
    specStruct.blackList = 1;
    
    switch lower(specStruct.dataset)
        case 'healthy1'
            specStruct.exerciseAcceptSuffixString = {'SLO'};
            specStruct.exerciseAcceptSuffixInstance = [0];
            
        otherwise
    end

    % load the file paths
    pathToRawDataAugmented = fullfile(pathToRawData, dataPackage{ind_package}.datasetSpecs.dataPathSuffix);
    fileStackTemp = loadPatientFilepaths(pathToRawDataAugmented, specStruct);
    
    for ind_fileStack = 1:length(fileStackTemp)
        currFilePath = fileStackTemp{ind_fileStack}.filePath;
        currExerciseName = fileStackTemp{ind_fileStack}.exerciseName;
        fileStackTemp{ind_fileStack}.datasetName = lower(specStruct.dataset);
        
        fileStackTemp{ind_fileStack}.exerciseCropLength = specStruct.exerciseCropLength;

        switch lower(specStruct.dataset)
            case {'healthy1', 'healthy2', 'stjoseph1', 'tri1'}
                ekfPath = '2015_03_23';
                fileStackTemp{ind_fileStack}.filepathDataSegCrop = fullfile(currFilePath, 'Segmentation_cropping', 'SegmentData_Cropping_Manual.header');
                fileStackTemp{ind_fileStack}.filepathDataSegManual = fullfile(currFilePath, specStruct.manSeg, 'SegmentData_Manual_Manual.header');

%         fileStackTemp{ind_fileStack}.headerMocap = fullfile(currFilePath, 'Cortex', 'MocapData_Cortex.header');
                fileStackTemp{ind_fileStack}.filepathDataImu{1} = fullfile(currFilePath, 'Shimmer', 'SensorData_LKnee.header');
                fileStackTemp{ind_fileStack}.filepathDataImu{2} = fullfile(currFilePath, 'Shimmer', 'SensorData_LAnkle.header');
                fileStackTemp{ind_fileStack}.filepathDataImu{3} = fullfile(currFilePath, 'Shimmer', 'SensorData_RKnee.header');
                fileStackTemp{ind_fileStack}.filepathDataImu{4} = fullfile(currFilePath, 'Shimmer', 'SensorData_RAnkle.header');
                fileStackTemp{ind_fileStack}.filepathDataEkf = fullfile(currFilePath, 'EKF', ekfPath, ['ekf.header']);
                
            case 'wojtusch2015_ichr'
                fileStackTemp{ind_fileStack}.filepathDataEkf = fullfile(currFilePath, 'EKF', 'EKF.mat');
                
            case 'kulic2009_tro'
                fileStackTemp{ind_fileStack}.filepathDataEkf = '';
                
            case 'myo'
                fileStackTemp{ind_fileStack}.filepathDataEkf = fullfile(currFilePath, 'Shapehand', 'kin.mat');
                
            case 'squats_tuat_2011'
                jointAngleFile = searchForFileByExt(fullfile(currFilePath, 'JointAngles', 'IK_2016-01_MK'), '*.mat');
                fileStackTemp{ind_fileStack}.filepathDataEkf = fullfile(currFilePath, 'JointAngles', 'IK_2016-01_MK', jointAngleFile);
                fileStackTemp{ind_fileStack}.filepathDataSegManual = fullfile(currFilePath, specStruct.manSeg, 'SegmentData_Manual_Manual.data');
                
            case 'squats_tuat_2015'
                fileStackTemp{ind_fileStack}.filepathDataSegCrop = fullfile(currFilePath, 'Segmentation_cropping', 'SegmentData_Cropping_Manual.header');
                                
                jointAngleFile = searchForFileByExt(currFilePath, 'Kinematics_meas*.mat');
                fileStackTemp{ind_fileStack}.filepathDataMocap = fullfile(currFilePath, jointAngleFile);
                fileStackTemp{ind_fileStack}.filepathDataEkf = fullfile(currFilePath, jointAngleFile);
                
                dynamicsFile = searchForFileByExt(currFilePath, 'Dynamics_meas_*.mat');
                fileStackTemp{ind_fileStack}.filepathDataGrf = fullfile(currFilePath, dynamicsFile);
                
                torqueFile = searchForFileByExt(currFilePath, 'Joint_torques*.mat');
                fileStackTemp{ind_fileStack}.filepathDataTorque = fullfile(currFilePath, torqueFile);
                
                fileStackTemp{ind_fileStack}.filepathDataSegManual = fullfile(currFilePath, specStruct.manSeg, 'SegmentData_Manual_Manual.data');
                
            case 'taiso_ut_2009'
                jointAngleFile = searchForFileByExt(fullfile(currFilePath, 'JointAngles', 'anm'), '*.anm');
                fileStackTemp{ind_fileStack}.filepathDataEkf = fullfile(currFilePath, 'JointAngles', 'anm', jointAngleFile);
                
                fileStackTemp{ind_fileStack}.filepathDataSegManual = fullfile(currFilePath, specStruct.manSeg, 'SegmentData_Manual_Manual.csv');
                  
%                 segFile = searchForFileByExt(fullfile(currFilePath, 'Segmentation_manual_GV'), 'Section*.m');
%                 fileStackTemp{ind_fileStack}.segmentGuidancePath = fullfile(fullfile(currFilePath, 'Segmentation_manual_GV'), segFile);
                
                dynamicsFile = searchForFileByExt(fullfile(currFilePath, 'Forces'), '*.forces');
                fileStackTemp{ind_fileStack}.filepathDataGrf = fullfile(currFilePath, 'Forces', dynamicsFile);
                
                
            case 'doppel'
                jointAngleFile = searchForFileByExt(currFilePath, '*doppel*.data');
                fileStackTemp{ind_fileStack}.filepathDataEkf = fullfile(currFilePath, jointAngleFile);
                
                segFile = searchForFileByExt(currFilePath, '*SegmentData_Manual*.data');
                fileStackTemp{ind_fileStack}.filepathDataSegManual = fullfile(currFilePath, segFile);
        end
         
%         fileStackTemp{ind_fileStack}.segmentPath = fullfile(currFilePath, specStruct.manSeg, 'SegmentData_Manual_Manual.data');
    end
   
    fileStack = [fileStack, fileStackTemp];
end

% iterate through the filestack and pull out all the exercise types we are
% including
exerciseType = [];
for ind_fileStack = 1:length(fileStack)
    if isempty(exerciseType)
        exerciseType{1} = fileStack{ind_fileStack}.exerciseType;
    elseif ~isempty(exerciseType) && sum(strcmpi(exerciseType, fileStack{ind_fileStack}.exerciseType)) == 0
        exerciseType{end+1} = fileStack{ind_fileStack}.exerciseType;
    end
end

for ind_fileStack = 1:length(exerciseType)
    switch exerciseType{ind_fileStack} % replace the taisou text with the actual contents of the taiso exercise
        case 'TAIS_STD_ONE'
            exerciseListFull = {...
                'MoveBody3Times', 'ArmCircleInnerOuter', 'SideBendingRightLeft', 'FrontBending', ...
                'WaistRotationRightLeft', 'LegsArms', 'TouchFootLeftRight', 'BigCirclesRightLeft', 'Jump'}; % TAISOU_1
            
            exerciseType{ind_fileStack} = exerciseListFull{1};
            exerciseType = [exerciseType exerciseListFull{2:end}];
    end
end

exerciseStruct.exerciseType = exerciseType;
exerciseStruct.count = length(exerciseType);
end