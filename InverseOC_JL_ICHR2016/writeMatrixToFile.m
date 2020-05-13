function writeMatrixToFile(fileID, outputString, matrixToWrite)
    
    matrixToWrite = matrixToWrite';
    fmt = cell(1, size(matrixToWrite, 1));
    fmt(:) = {',%f'};
    matrixOutputStr = horzcat([',,' fmt{:} '\n']);

    fprintf(fileID, outputString);
    fprintf(fileID, matrixOutputStr, matrixToWrite);
