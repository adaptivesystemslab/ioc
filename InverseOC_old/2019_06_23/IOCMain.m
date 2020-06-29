%clear all;
clc;

%% Add paths to directories with model definition and util functions
addpath(genpath('../Utils/')); 
addpath(genpath('../Logic/'));
addpath(genpath('../Libraries/rvctools'));
addpath(genpath('../Libraries/cvx'));
addpath(genpath('../Libraries/rl/rlMatlabWrapper'));
addpath(genpath('../Libraries/rl/Common')); 
       
%% Add cvx to path
if ~exist('cvx_begin', 'file')
    cvx_setup;
end

%% Add Peter Corke's toolbox to library path if required
if ~exist('rtbdemo', 'file')
    startup_rvc;
end

%% Load json file with list of all trials on which IOC will be run
status = 1;

%% Create and/or look for folder where solutions are going to be saved
currentDate = erase(string(datetime("today")),"-");
savePath = sprintf("../Data/IOC/%s/", currentDate);

if ~exist(savePath, 'dir')
   status = mkdir(char(savePath)); 
end

if status
    %% Load json with information about each trial
    configFile = jsondecode(fileread('../Data/IOC_CostAnalysis.json'));
    
    n = length(configFile.Files);
    
    for i=1:n
        trial = configFile.Files(i);
        fprintf("Processing %s file \n", trial.name);
        runIOC(trial, savePath, 1);
    end
else
    disp("Error finding/creating path where results will be stored \n Nothing to do here");
end