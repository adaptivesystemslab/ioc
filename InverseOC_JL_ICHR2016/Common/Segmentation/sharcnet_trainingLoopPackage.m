function package = ...
    sharcnet_trainingLoopPackage(dataset, patientList, ind_patient, ...
    exerciseList, ind_exercise, param, exerciseCropLength)

    package.dataset = dataset;
    package.patient = patientList{ind_patient};
    try
        package.session = exerciseList.sessionList{ind_exercise};
    catch
        package.session = [];
    end
    
    package.exerciseAcceptPrefixString = exerciseList.exerciseString{ind_exercise};
    package.exerciseAcceptPrefixInstance = exerciseList.exerciseInstance{ind_exercise};

    package.datasetSpecs = datasetSpecs(dataset);
    package.exerciseCropLength = exerciseCropLength; % firsthalf, secondhalf, full
    
    package = sharcnet_paramSet(package, param);
end