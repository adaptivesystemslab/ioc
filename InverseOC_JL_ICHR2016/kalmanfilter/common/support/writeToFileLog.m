function writeToFileLog(basepathTarget, currFileEntry, messagePrefix, messageSuffix, algorithmParam)
    logPath = globalConstants_filepaths.fileFullFileRuntimeLog(basepathTarget);
    checkMkdir(logPath);
    
    [fId, err] = fopen(logPath, 'a');
    if ~isempty(currFileEntry) && isfield(currFileEntry, 'exerciseName')
        fprintf(fId, '%s: Subject %s, Exercise %s - %s\n', ...
            messagePrefix, currFileEntry.subjectString, currFileEntry.exerciseName, ...
            messageSuffix);
    elseif ~isempty(currFileEntry)
        fprintf(fId, '%s: Subject %s - %s\n', ...
            messagePrefix, currFileEntry.subjectString, ...
            messageSuffix);
    else
        fprintf(fId, '%s: %s\n', ...
            messagePrefix, ...
            messageSuffix);
    end

    fclose(fId);
end 