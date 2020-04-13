function csv_populate(matData, masterPathCsv)
    if ~exist(masterPathCsv, 'file')
        %             if doesn't exist write header
        header = 'Subject,Model,CF,Direction,DataFile';
        header = [header ',timeElapsedTotal,numFrames,timeElapsedPerLength'];
        header = [header '\n'];
    else
        header = '';
    end

    [fileID, errmsg] = fopen(masterPathCsv, 'a');
    fprintf(fileID, header);

    % write data to log file
    runNameSplit = strsplit(matData.trialInfo.runName, '_');
    subpath = strsplit(matData.trialInfo.subpath, '/');
    fprintf(fileID,'%s,%s,%s,%s,%s',runNameSplit{1}, runNameSplit{2}, ...
        runNameSplit{3}, matData.trialInfo.hWinAdvFlag, subpath{end});
    
%     framesProc = (matData.frameInds(end) - matData.frameInds(1) + 1);
    framesProc = length(matData.progress);
    timeElapsedPerFrame = matData.timeElapsed / framesProc;
    fprintf(fileID,',%f,%f,%f',matData.timeElapsed, framesProc, timeElapsedPerFrame);
    
    fprintf(fileID,'\n');
    fclose(fileID);
end