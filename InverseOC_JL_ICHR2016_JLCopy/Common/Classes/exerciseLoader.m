function exerciseStruct = exerciseLoader(varargin)
    % exerciseLoader loads all associated data for a given dataset, given
    % user specifications. This function would need to be updated to each
    % local machine to ensure that the pathing is correct. For an example,
    % please see exerciseLoader_mwe
    
    specStruct = varargin{1};
    pathToRawData = callDatasetBasePath(specStruct); % update the contents of this function with local specifications

    % iterate and load all patient/session/exercise combination that meets
    % user requirements
    fileStack = loadPatientFilepaths(pathToRawData, specStruct);

    % load data
    exerciseStruct = cell(size(fileStack));
    for ind_fileStack = 1:length(fileStack)
        exerciseStruct{ind_fileStack} = exerciseDataHandle(fileStack{ind_fileStack}.filePath, varargin{2:end});
        exerciseStruct{ind_fileStack}.load;
    end
end