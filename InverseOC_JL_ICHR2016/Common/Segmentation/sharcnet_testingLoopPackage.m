function package = ...
    sharcnet_testingLoopPackage(dataset, patientList_training, patientList_full, ...
     exerciseList, ind_exercise, param, exerciseCropLength, specialCommand, specialPackage)

     if iscell(patientList_full)
         % the testing array may be passed in as a cell array. if that's
         % the case, convert it from the cell array to a normal array first
         patientList_full = patientList_full{1};
     end
 
    package.dataset = dataset;
    package.patient = setxor(patientList_full, patientList_training.patient);
    
    % there are no intersections
    if length(package.patient) > length(patientList_full) || isempty(package.patient)
        % in that case, just keep the original entries as the patient to
        % load
        package.patient = patientList_full;
    end
    
    % session for that patient
    for i = 1:length(exerciseList.sessionList)
        if ~isempty(exerciseList.sessionList{i})
            package.session{i} = setxor(patientList_training.session{i}, exerciseList.sessionList{i});
        end
    end
    
    package.exerciseAcceptPrefixString = exerciseList.exerciseString{ind_exercise};
    package.exerciseAcceptPrefixInstance = exerciseList.exerciseInstance{ind_exercise};
    
    package.datasetSpecs = datasetSpecs(dataset);
    package.exerciseCropLength = exerciseCropLength; % firsthalf, secondhalf, full
    
    package = sharcnet_paramSet(package, param);    
    
    if exist('specialCommand', 'var')
        switch specialCommand
            case 'exercisesetxor'
                % if a specific instance number is included, then attach it
                % to the prefix string
                acceptPrefixString_package = package.exerciseAcceptPrefixString;
                acceptPrefixInstance_package = package.exerciseAcceptPrefixInstance;
                
                acceptPrefixString_specialPackage = specialPackage.exerciseAcceptPrefixString;
                acceptPrefixInstance_specialPackage = specialPackage.exerciseAcceptPrefixInstance;
                
                for ind = 1:length(acceptPrefixString_package)
                    if acceptPrefixInstance_package(ind) > 0
                        acceptPrefixString_package{ind} = [acceptPrefixString_package{ind} num2str(acceptPrefixInstance_package(ind))];
                    end
                end
                
                for ind = 1:length(acceptPrefixString_specialPackage)
                    if acceptPrefixInstance_specialPackage(ind) > 0
                        acceptPrefixString_specialPackage{ind} = [acceptPrefixString_specialPackage{ind} num2str(acceptPrefixInstance_specialPackage(ind))];
                    end
                end
                
                [acceptprefixstring, a, b] = setxor(acceptPrefixString_package, acceptPrefixString_specialPackage);
                package.exerciseAcceptPrefixString = package.exerciseAcceptPrefixString(a);
                package.exerciseAcceptPrefixInstance = package.exerciseAcceptPrefixInstance(a);
        end
    end
end

