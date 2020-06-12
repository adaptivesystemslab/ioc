function [fileName] = searchForFileByExt(filePathSearch, fileNameSearch)
% searches filePath for all files of a given extension, and return the
% filename. if multiple files are found, the first one is returned. 

pathToSearch = fullfile(filePathSearch, fileNameSearch);
filePathDir = dir(pathToSearch);

if ~isempty(filePathDir)
    if length(filePathDir) == 1
         fileName = filePathDir.name;
    else
         fileName = filePathDir{1}.name;
    end
else
    fileName = [];
end