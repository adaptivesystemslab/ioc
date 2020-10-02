function [checkIk, checkFk] = checkFSOS(currFileEntry, filepathInstance, algorithmParam)
    pathToSave_featureSetEkf_ik = globalConstants_filepaths.fileFullEkfIk(filepathInstance, currFileEntry, algorithmParam);
    pathToSave_featureSetEkf_fk = globalConstants_filepaths.fileFullEkfFk(filepathInstance, currFileEntry, algorithmParam);

    if ~exist(pathToSave_featureSetEkf_ik, 'file') && ~exist(pathToSave_featureSetEkf_fk, 'file')
        checkIk = 0;
        checkFk = 0;
    elseif exist(pathToSave_featureSetEkf_ik, 'file')
        checkIk = 1;
        checkFk = 0;
    else
        checkIk = 0;
        checkFk = 1;
    end
end