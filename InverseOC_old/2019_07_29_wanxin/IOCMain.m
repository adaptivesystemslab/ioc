% test RL-based cost function derivations with the original files generated
% by PC

%clear all;
clc;
tic;

%% Add paths to directories with model definition and util functions
setPaths();

%% Define internal settings
runParam.displayInfo = 'final'; % 'none', 'final', 'verbose'
runParam.saveIntermediate = 200;
runParam.gamma = 10;
runParam.dimWeights = 5;
runParam.maxWinLen = 100;
runParam.start = 0;
% runParam.configFile = '../Data/IOC_testCostFunction.json';
runParam.configFile = '../Data/IOC_IITFatigue_WJTest.json';
% runParam.configFile = '../Data/IOC_IITFatigue_full.json';
% runParam.configFile = '../Data/IOC_IITFatigue_full_3dof.json';

% these variables are not set, but saved for record keeping purposes
runParam.windowWidth = IOCInstanceNew.winSize;
runParam.delta = IOCInstanceNew.delta;

%% Load json file with list of all trials on which IOC will be run
configFile = jsondecode(fileread(runParam.configFile));

%% Create and/or look for folder where solutions are going to be saved
currentDate = datestr(datetime("now"),"yyyy_mm_dd_HH_MM_SS");
savePath = sprintf("../Data/IOC/%s/", currentDate);

if checkMkdir(char(savePath))
    %% Load json with information about each trial
    n = length(configFile.Files);
    
    for i=1:n
        if length(configFile.Files) == 1
            trial = configFile.Files;
        else
            trial = configFile.Files(i);
            
            if iscell(trial)
                trial = trial{1};
            end
        end
        
        % if the source matfile is not found in the json path, search these
        % following locations as well
        potentialBasePaths = {'/project/6001934/data/', ...
            trial.basepath};
        
        for j = 1:length(potentialBasePaths)
            targetPath = fullfile(potentialBasePaths{j}, trial.subpath);
            
            if exist(targetPath, 'file')
                break;
            end
        end
        
        fprintf("Processing %s file \n", trial.name);
        trial.path = targetPath
        
        currentDate = datestr(datetime("now"),"yyyy_mm_dd_HH_MM_SS");
        savePath = sprintf("../Data/IOC/%s/", currentDate);
        checkMkdir(savePath);
        IOCRun_WJTest2(trial, savePath, runParam);
        
%         currentDate = datestr(datetime("now"),"yyyy_mm_dd_HH_MM_SS");
%         savePath = sprintf("../Data/IOC/%s/", currentDate);
%         checkMkdir(savePath);
%         IOCRun_20190724(trial, savePath, runParam);
%         
%         currentDate = datestr(datetime("now"),"yyyy_mm_dd_HH_MM_SS");
%         savePath = sprintf("../Data/IOC/%s/", currentDate);
%         checkMkdir(savePath);
%         IOCRun(trial, savePath, runParam);
    end
else
    disp("Error finding/creating path where results will be stored \n Nothing to do here");
end

toc