function patientStruct = patientLoader(varargin)
    % exerciseLoader loads all associated data for a given dataset, given
    % user specifications. This function would need to be updated to each
    % local machine to ensure that the pathing is correct. For an example,
    % please see exerciseLoader_mwe
    
    specStruct = varargin{1};
    options = varargin{2};
    
    pathToRawData = callDatasetBasePath(specStruct);

    if isempty(specStruct.patient)
        % iterate though the dir and pull out all the patients
        fprintf('Please populate specStruct.patient a patient or a range of patients. \n');
        patientStruct = [];
        return
    end
    
    patientArray = specStruct.patient;
    patientStruct = patientDataHandle.empty(numel(patientArray),0);
    for ind_patient = 1:length(patientArray)
        specStruct.patient = patientArray(ind_patient);
        patientStruct(ind_patient) = patientDataHandle(pathToRawData, specStruct, options);
    end
end