function normArray = calc_normalization_rangeArray
    
    normArray = [];

    for i = 2:7
        outputString = ['win2011p' num2str(i) '_3r_10_61_20_startmid1endonknot_8_none_100_none_rang']; % win_dof_knot_win_shift_constraint
        [specStruct, variableFactorsInput] = script2_testingSpec(outputString, 0);
        normArrayNew = calcNormalization(specStruct, variableFactorsInput);
        
        normArray = [normArray; normArrayNew];
    end

    normMean = mean(normArray, 1);
end

function normArray = calcNormalization(specStruct, variableFactorsInput)    
    specStruct.datasetSpecs = datasetSpecs(specStruct.dataset);
    pathToRawData = fullfile('C:', 'Documents', 'aslab', 'data', specStruct.datasetSpecs.dataPathSuffix);
    runSettings.variableFactors = variableFactorsInput;
    
    fileStackTemp = loadPatientFilepaths(pathToRawData, specStruct);
    
    for ind_fileStack = 1:length(fileStackTemp)
        currFilestack = fileStackTemp{ind_fileStack};
        currFilePath = currFilestack.filePath;
        currFilestack.dataset = lower(specStruct.dataset);
        
        switch currFilestack.dataset
            case 'squats_tuat_2011'
                manSeg = 'Segmentation_manual_MK';
                
                jointAngleFile = searchForFileByExt(fullfile(currFilePath, 'JointAngles', 'IK_2016-10_JL'), 'jointAngles*.mat');
                filesToLoad{1} = fullfile(fullfile(currFilePath, 'JointAngles', 'IK_2016-10_JL'), jointAngleFile);
                
                mocapFile = searchForFileByExt(fullfile(currFilePath, '..', 'prm'), 'sub*.mat');
                filesToLoad{2} = fullfile(fullfile(currFilePath, '..', 'prm'), mocapFile);
                
                filesToLoad{3} = 'C:\Documents\aslab\data\Squats_TUAT_2015-12\Subject01\Session1\SQUA_STD_NON1\Dynamics_meas_Sb1_Tr1.mat';
                
            case 'squats_tuat_2015'
                manSeg = 'Segmentation_manual_JL';
                
                jointAngleFile = searchForFileByExt(currFilePath, 'Kinematics_meas*.mat');
                filesToLoad{1} = fullfile(currFilePath, jointAngleFile);
                
                dynamicsFile = searchForFileByExt(currFilePath, 'Dynamics_meas_*.mat');
                filesToLoad{2} = fullfile(currFilePath, dynamicsFile);
        end
        
        manSegLoadPath = fullfile(currFilePath, manSeg, 'SegmentData_Manual_Manual.data');
        
        % load the data to know how much of the file needs parsing
        param = setup_main(filesToLoad, manSegLoadPath, currFilestack, 'win', runSettings, [], []);
        %             param.coeff_cf.array = [param.coeff_cf.ddq param.coeff_cf.dddq param.coeff_cf.ddx param.coeff_cf.dddx param.coeff_cf.tau param.coeff_cf.dtau param.coeff_cf.ddtau param.coeff_cf.ep param.coeff_cf.ek param.coeff_cf.geo param.coeff_cf.en param.coeff_cf.cop param.coeff_cf.dcop];
       
        normArray = param.coeff_cf.array;
    end
end