function fId = fopenCheck(exportPath, header)
    fId = -1;
    
    checkCounter = 0;
    
    while fId < 0
        fId = fopen(exportPath, 'a');
        
        if fId < 0
            pause(3); % couldn't open it because something else has access to it
            checkCounter = checkCounter + 1;
            
            if checkCounter > 500
                error('fopenCheck: %s cannot be opened', exportPath);
            end
        end
    end
    
    if ~isempty(header)
        fprintf(fId, '%s\n', header);
    end
end