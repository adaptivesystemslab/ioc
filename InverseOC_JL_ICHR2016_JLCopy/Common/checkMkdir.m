function checkMkdir(targetFolder)
    if ~exist(targetFolder, 'dir')
        mkdir(targetFolder);
    end
end