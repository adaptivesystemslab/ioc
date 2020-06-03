function [currTrainingPxStr, currTestingPxStr, pxStr] = setupOutputConfigNames(currTrainingPackage, currTestingPackage, dataSelect)

currTrainingPxStr = '';
for ind_str = 1:length(currTrainingPackage{1}.patient)
    currTrainingPxStr = [currTrainingPxStr '_' num2str(currTrainingPackage{1}.patient(ind_str))];
end

if length(currTrainingPxStr) > 10
    currTrainingPxStr = currTrainingPxStr(1:10);
end

% set up the training/testing packages
if ~isempty(currTestingPackage)
    currTestingPxStr = '';
    for ind_str = 1:length(currTestingPackage{1}.patient)
        currTestingPxStr = [currTestingPxStr '_' num2str(currTestingPackage{1}.patient(ind_str))];
    end
    
    if length(currTestingPxStr) > 10
        currTestingPxStr = currTestingPxStr(1:10);
    end
    
    pxStr = ['tr' currTrainingPxStr ...
        '-te' currTestingPackage{1}.dataset(1:3) '-' currTestingPxStr ...
        '-co' dataSelect.settings.settingName];
else
    currTestingPxStr = 'NoTestingData';
    pxStr = 'NoTestingData';
end

if length(pxStr) > 25
    pxStr = pxStr(1:25);
end