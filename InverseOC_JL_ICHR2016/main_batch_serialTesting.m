% iterative version
knotTest = 1:30;
for ind_knotTest = 1:length(knotTest)
pathToRawData = 'C:\Documents\aslab\data\Squats_TUAT_2015-12\';
specStruct.dataset = 'Squats_TUAT_2015';
specStruct.patient = 1; %[1 2 3 4 5 6 7 8];
specStruct.exerciseAcceptPrefixString = {'SQUA_STD'};
specStruct.exerciseAcceptPrefixInstance = 1;
specStruct.datasetSpecs = datasetSpecs(specStruct.dataset);
specStruct.manSeg = 'Segmentation_manual_JL';

fileStackTemp = loadPatientFilepaths(pathToRawData, specStruct);

outputBasePath = 'C:\Documents\MATLABResults\IOCProject\realtesting\';
outputInstancePath = datestr(now, 'yyyy-mm-dd-HH-MM-SS');
outputPath = [outputBasePath outputInstancePath];
checkMkdir(outputPath);

for ind_fileStack = 1:length(fileStackTemp)
    currFilestack = fileStackTemp{ind_fileStack};
    currFilePath = currFilestack.filePath;
    currInstName = ['Subj' num2str(currFilestack.subjectNumber) '_' currFilestack.exerciseName];        
    
    jointAngleFile = searchForFileByExt(currFilePath, 'Kinematics_meas*.mat');
    trajToLoadPath = fullfile(currFilePath, jointAngleFile);
    
    dynamicsFile = searchForFileByExt(currFilePath, 'Dynamics_meas_*.mat');
    dynToLoadPath = fullfile(currFilePath, dynamicsFile);
    
    manSegLoadPath = fullfile(currFilePath, specStruct.manSeg, 'SegmentData_Manual_Manual.data');

    main(trajToLoadPath, dynToLoadPath, manSegLoadPath, outputPath, currInstName, knotTest(ind_knotTest));
end
end
