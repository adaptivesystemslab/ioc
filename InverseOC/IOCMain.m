% test RL-based cost function derivations with the original files generated by PC

clear all;
close all;
clc;
tic;

%% Define internal settings
% Jonathan's Configurations
% configFilePath = '../Data/IOC_gitupload_test.json';
% configFilePath = '../Data_json/CarrenoConfig/IOC_ExpressiveData_test.json';
% configFilePath = '../Data_json/LinConfig/IOC_IITFatigue_test.json';
% configFilePath = '../Data_json/LinConfig/IOC_Healthy1.json';
% configFilePath = '../Data_json/LinConfig/IOC_Jumping2D.json';
% configFilePath = '../Data_json/IOC_IITFatigue_test.json';

% Wanxin's Configurations
configFilePath = '../Data_json/JinConfig/Squat_IIT/IOC_IITFatigue_Test_Sub1.json';
savePath='../Data/JinIOCResults/IOC_Squat_IIT/Sub1/';
% configFilePath = '../Data_json/JinConfig/KneeHipFlexion/IOC_KHFHealthy_Sub1.json';
% savePath='../Data/JinIOCResults/IOC_KneeHipFlexion/Sub1/';
% configFilePath = '../Data_json/JinConfig/HipFlexion/IOC_HFHealthy_Sub1.json';
% configFilePath = '../Data_json/JinConfig/Squat/IOC_SQUHealthy_Sub1.json';
% configFilePath = '../Data_json/JinConfig/Sit/IOC_SITHealthy_Sub1.json';
%% Create and/or look for folder where solutions are going to be saved
overwriteFiles = 1;

%% Set up internal parameters
% Add paths to directories with model definition and util functions
setPaths();

% Load json file with list of all trials on which IOC will be run
configFile = jsondecode(fileread(configFilePath));
n = length(configFile.Files);

for i=3
    runParam = [];
    configFileI = configFile.Files(i);

    % if the source matfile is not found in the json path, search these
    % following locations as well, such as for Sharcnet deployment
    potentialBasePaths = {'/project/6001934/data/', ...
        configFileI.basepath, ...
        'H:/data'};
    
    % load the specific trialinfo
    fprintf("Processing %s file \n", configFileI.runName);
    trialInfo = loadTrialInfo(configFileI, configFile, potentialBasePaths, configFilePath);
    
    % does the target folder already exist? 
    [status, alreadyExist] = checkMkdir(savePath);
    
    if ~alreadyExist || overwriteFiles
%         IOCRun(trialInfo, subsavePath);
        IOCIncomplete(trialInfo,savePath)
    end
end

toc