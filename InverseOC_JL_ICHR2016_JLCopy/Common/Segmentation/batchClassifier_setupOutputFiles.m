[currTrainingPxStr, currTestingPxStr, pxStr] = setupOutputConfigNames(currTrainingPackage, currTestingPackage, dataSelect);

exportSuffix = [dataSelect.settings.settingName '_' dimReductSelect.settingName '_' classifierSelect.settingName '_' aggregatorSelect.settingName];

if length(exportSuffix) > 25 % this prevents the filepath strings from getting too long, which can cause MATLAB to error out
    exportSuffix = exportSuffix(1:25);
end

% set up the actual file paths for the various outputs
exportPrefix = fullfile(dataSelect.exportBasePath, dataSelect.outputBase, batchInstancePath, batchSettings.exportPathSuffix);
exportPath = fullfile(exportPrefix, pxStr);

overallSummaryTestingPath  = fullfile(exportPrefix, [batchInstancePath '_summary_overall.csv']);
instanceSummaryTestingPath = fullfile(exportPrefix, [batchInstancePath '_summary_instance.csv']);
overallCSVTrainingPath     = fullfile(exportPrefix, ['score_' pxStr '_Training.csv']);
overallCSVTestingPath      = fullfile(exportPrefix, ['score_' pxStr '_Testing.csv']);
exportPathInstance         = fullfile(exportPath, exportSuffix);
instanceDiary              = [exportPathInstance '_diary.txt'];
instanceCSVTrainingPath    = [exportPathInstance '_Training.csv'];
instanceCSVTestingPath     = [exportPathInstance '_Testing.csv'];

checkMkdir(exportPath);
checkMkdir(exportPathInstance);