function createGlobalHeader(targetFile)
    % make header
    fId = fopen(targetFile, 'w');
    
    fprintf(fId, '<?xml version="1.0"?>\n');
    
    fclose(fId);
end