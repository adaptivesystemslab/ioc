function [status, alreadyExist] = checkMkdir(targetFolder)
    % check to see the target folder exist. if not, create the folder
    
    % if the target path passed in is actually a file path, strip out the file 
    [targetFolder2,name,ext] = fileparts(targetFolder);
    if ~isempty(ext)
        targetFolder = targetFolder2;
    end
    
    alreadyExist = exist(targetFolder, 'dir');
    status = 1;
    if ~alreadyExist
        [status,msg] = mkdir(targetFolder);
    end
end