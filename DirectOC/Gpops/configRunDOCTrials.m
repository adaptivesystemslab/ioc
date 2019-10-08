%% Load file with configuration for different trials and run one by one

clc;clear;close all;

if ~exist('gpops2', 'file')
    cd ../../gpops/GPOPS-II/;
    gpopsMatlabPathSetup;
    cd ../../expressive-ioc/DirectOC/;
end

% Add paths to directories with model definition and util functions
addpath(genpath('../Utils/')); 
addpath(genpath('../Logic/'));
addpath(genpath('../InverseOC/'));
addpath(genpath('../Analysis/'));
addpath(genpath('../Libraries/rl/rlMatlabWrapper'));
addpath(genpath('../Libraries/rl/Common'));
addpath(genpath('../../gpops')); 

% Define auxiliary variables
figurePath = "";
N = 4;
status = 1;
plot = 0;

% Create and/or look for folder where solutions are going to be saved
currentDate = erase(string(datetime("today")),"-");
savePath = sprintf("../Data/OC/%s/", currentDate);

if ~exist(savePath, 'dir')
   status = mkdir(char(savePath)); 
end


% Generate labels for plots
controlNames = [""]; angleNames = [""]; velNames = [""];
for i=1:N
        controlNames(:,end+1) = sprintf("u_%d", i);
        angleNames(:,end+1) = sprintf("\\theta_%d", i);
        velNames(:,end+1) = sprintf("\\dot\\theta_%d", i);
end

% Create structure that will summarize trial information necessary for IOC
% This structure will be parsed and saved as a json file
jsonStruct = struct;

if status
    % Load json with information about each trial
    configFile = jsondecode(fileread('../Data/DOC_ListFiles_NewFormat.json'));
    numTrials = length(configFile.Trials);
    
    for i=3:3
        
        % Create setup and define continuous functions for GPOPS2
        currentTrial = configFile.Trials(i);
        weights = currentTrial.weights;
                
        % Create dynamic model and IOC instance associated to it. This
        % instance will be used to compute target features
        dynModel = getModel(currentTrial.model, currentTrial.modelType);
        ioc = IOCInstanceNew(dynModel, 0.01);
        ioc.init(currentTrial);
                
        % Define optimization setup based on json entry
        % robotContinuous = @(input)parameterizedDynamics(input, weights, dynModel);
        robotContinuous = @(input)parameterizedDynamics(input, weights, ioc);
        
        % Check type of initial guess: manual definition or spline
        % computation
        if ischar(currentTrial.guess)
           [q, dq, tau, time] = getSplineGuess(dynModel, currentTrial.initialState,...
               currentTrial.finalState);
           
           tempStruct = struct('jointAngles', q, 'angularVelocities', dq,...
               'control', tau, 'time', time);
           currentTrial.guess = tempStruct;
        end
        
        setup = defineGPOPSOptimization(currentTrial, robotContinuous); 
        
        %  Solve problem using GPOPS2
        output = gpops2(setup);

        % Extract Solution and save it in data folder
        solution = output.result.solution;
        time = solution.phase.time;
        state = solution.phase.state;
        control = solution.phase.control;
        mat = [time, control, state];
        finalPath = sprintf("%s%s_%s.mat", savePath, currentTrial.name,...
             currentTrial.model);
        save(char(finalPath), 'mat');
        
        % Add information to json file structure
        jsonField.name = string(currentTrial.name);
        jsonField.path = finalPath;
        jsonField.weights = weights;
        jsonField.candidateFeatures = ["TorqueSquared", ""];
        jsonField.model = string(currentTrial.model);
        jsonField.modelType = string(currentTrial.modelType);
        jsonStruct.Files(i) = jsonField;        
        
        % Plot Solution.
        if plot
            drawSavePlots(time', state(:,1:size(control,2)), angleNames(:, 2:N+1), "", "");
            drawSavePlots(time', state(:,size(control,2)+1:end), velNames(:, 2:N+1), "", "");
            drawSavePlots(time', control, controlNames(:, 2:N+1), "", "");
            pause;
            close all;
        end

    end
    
    % Parse, format and save json file
    json = jsonencode(jsonStruct);
    json = strrep(json, ',', sprintf(',\r'));
    json = strrep(json, '[{', sprintf('[\r{\r'));
    json = strrep(json, '}]', sprintf('\r}\r]'));
    
    fid = fopen('../Data/IOC_FilesList.json', 'w');
    if fid == -1
        error('Cannot create JSON file');
    end
    fwrite(fid, json, 'char');
    fclose(fid);
    
else
    fprintf("Error on identification of data's destination folder");
end


