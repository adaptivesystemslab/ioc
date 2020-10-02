function calc_calibration_matrices(ind_subj, databaseSpecTable, basepathSource, subfolderMocap, subfolderImu, filepathCalibrationTarget, fieldsImuModalityCalibration)
    % TODO: use currSubjId and dataSpecTable to load the proper IMU files, then
    % iterate though all the calibration techniques and determine which ones
    % need an external set of calibration matrices calculated. 
    
    for i = 1:length(fieldsImuModalityCalibration)
        currCalibration = fieldsImuModalityCalibration{i};
        currSubjId = databaseSpecTable.SubjectName{ind_subj};
        
        % build the save paths
        pathToSaveSuffix_calibStruct = ['imuCalib' '_' currSubjId '_'  currCalibration];
        savePath = fullfile(filepathCalibrationTarget, pathToSaveSuffix_calibStruct); % save the calibration structure to this
        
        imuTransforms = [];
        switch currCalibration
            case 'luinge2008'
                currExerciseName = 'HEFO_STD_NON1';
                algorithmParam = [];
                currFileEntry = setupCurrFileEntry_ssa(basepathSource, currExerciseName, subfolderMocap, subfolderImu, ind_subj, currSubjId, databaseSpecTable);
                
                modelInstance = rlModelInstance_ssa(ind_subj);
                dataInstance_hefo = rlDataInstance_imu(modelInstance);
                dataInstance_hefo.loadData(currFileEntry, algorithmParam);
                
                transform = luinge2008calib(dataInstance_hefo, modelInstance);
                
            case 'galinski2012'
                % TODO
                transform = galinski2012evaluation(dataInstance_hefo, modelInstance);
                
            case 'kong2016'
                % TODO
                transform = kong2016anatomical(dataInstance_hefo, modelInstance);
                
            case 'joukov2017'
                % TODO
                
            case 'seel2008'
                % TODO
                
            case 'chardonnens2012'
                % TODO
                transform = chardonnens2012effortless(dataInstance_hefo, modelInstance);
                
            case 'deVries2009'
                transform = devries2009magnetic(dataInstance_hefo, modelInstance);
        end
        
        save(savePath, 'imuTransforms');
    end
end