function [featureSet, pathToSave_featureSetEkf_ik, pathToSave_featureSetEkf_fk] = loadFSOS(currFileEntry, filepathInstance, algorithmParam)
    pathToSave_featureSetEkf_ik = globalConstants_filepaths.fileFullEkfIk(filepathInstance, currFileEntry, algorithmParam);
    pathToSave_featureSetEkf_fk = globalConstants_filepaths.fileFullEkfFk(filepathInstance, currFileEntry, algorithmParam);
    
    if ~exist(pathToSave_featureSetEkf_ik, 'file') && ~exist(pathToSave_featureSetEkf_fk, 'file')
        featureSet = [];
    elseif exist(pathToSave_featureSetEkf_ik, 'file')
        featureSet = rlFeatureSet_ioc();
        
        try
            featureSet.loadDataFromFile(pathToSave_featureSetEkf_ik);
        catch err
            featureSet.loadDataFromFileJointsOnly(pathToSave_featureSetEkf_ik);
        end
    else
        featureSet = rlFeatureSet_ioc();
        
        try
            featureSet.loadDataFromFile(pathToSave_featureSetEkf_fk);
        catch err
            featureSet.loadDataFromFileJointsOnly(pathToSave_featureSetEkf_fk);
        end
    end
end