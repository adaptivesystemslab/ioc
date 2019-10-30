%% First version of DOC code for RL model and Recovery matrix IOC

clear all;
clc;
tic;

%% Set up internal parameters
% Add paths to directories with model definition and util functions
setPaths();

potentialDataPaths = {'/project/6001934/data/', 'H:/data', '../../motionData/'...,
    'D:/aslab/data_IK'};

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                Read trial data and weights                              %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

status = 1;

% Create and/or look for folder where solutions are going to be saved
currentDate = erase(string(datetime("today")),"-");
savePath = sprintf("../../Data/OC/%s/", currentDate);

if ~exist(savePath, 'dir')
   status = mkdir(char(savePath)); 
end

frameInitPose = 1;
frameFinalPose = 100;

timeInit = 0;
timeFinal = 1;

if status
    % Load data and weights with hardcoded paths for now
    motion = load('/home/pcarreno/Documents/WaterlooPostDoc/Code/RecoveryMatrixIOC/outputData/Subj1_3DOF_3CF/data_000001_000999.mat');
    iocOutput = load('/home/pcarreno/Documents/WaterlooPostDoc/Code/RecoveryMatrixIOC/outputData/Subj1_3DOF_3CF/weights_000001_000999.mat');
    weightTrajectory = load('/home/pcarreno/Documents/WaterlooPostDoc/Code/RecoveryMatrixIOC/outputData/Subj1_3DOF_3CF/mat_results_cumulativeAllPass_result01_Subj1_3DOF_3CF_3DOF_3CF_Forward');
    
    motion = motion.outputVar_data;
    iocOutput = iocOutput.outputVar_weights;
    weightTrajectory = weightTrajectory.matSave;
    
    % Get information about model and features
    trialInfo = iocOutput.trialInfo;
    
    % Reset path to mat file if needed
    for j = 1:length(potentialDataPaths)
        targetPath = fullfile(potentialDataPaths{j}, trialInfo.subpath);
        if exist(targetPath, 'file')
            break;
        end
    end
    trialInfo.path = targetPath;
        
    % Create model and instantiate ioc instance
    model = getModel(trialInfo);
    dt = motion.dt;
    ioc = IOCInstance(model, dt);
    ioc.init(trialInfo);

    % Get bounds on state and control from motion file
    [stateBounds, controlBounds] = getBounds(motion);
    
    % Get initial and final state. Generate initial guess from spline
    initialState = [motion.q(frameInitPose,:), motion.dq(frameInitPose,:)];
    finalState = [motion.q(frameFinalPose,:), motion.dq(frameFinalPose,:)];
    
    guess = guessFromSpline(model, initialState, finalState);
    
    % Get weights to be used during objective computation. 
    % For now I assume weights correspond to avg recovery at finalState frame
    weights = weightTrajectory.weights(frameFinalPose,:);
    
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                      Set up function handles                            %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,x,u)( Dynamics(x,u,ioc) );
problem.func.pathObj = @(t,x,u)( Objective(t,x,u,weights,ioc) );


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                 Set up bounds on state and control                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.bounds.initialTime.low = timeInit;
problem.bounds.initialTime.upp = timeInit;
problem.bounds.finalTime.low = timeFinal;
problem.bounds.finalTime.upp = timeFinal;

problem.bounds.state.low = stateBounds.low;
problem.bounds.state.upp = stateBounds.upp;

problem.bounds.initialState.low = initialState';
problem.bounds.initialState.upp = initialState';

problem.bounds.finalState.low = finalState';
problem.bounds.finalState.upp = finalState';

problem.bounds.control.low = controlBounds.low;
problem.bounds.control.upp = controlBounds.upp;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                 Initialize trajectory with guess                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
problem.guess.time = guess.time;
problem.guess.state = [guess.q, guess.dq]';
problem.guess.control = guess.tau';


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                      Options for Transcription                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.options(1).nlpOpt = optimset('display','iter',...
    'MaxFunEval',1e5,...
    'tolFun',1e-6);
problem.options(1).method = 'rungeKutta';
problem.options(1).defaultAccuracy = 'low';

%problem.options.method = 'rungeKutta';
%problem.options.defaultAccuracy = 'low';
        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Solve!                                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

soln = optimTraj(problem);
toc
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Display the solution                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

figure(1); clf; hold on;

% t = linspace(soln.grid.time(1), soln.grid.time(end), 150);
% z = soln.interp.state(t);
% x = z(1,:);
% y = z(2,:);
% th = z(3,:);
% u = soln.interp.control(t);
% 
% tGrid = soln.grid.time;
% xGrid = soln.grid.state(1,:);
% yGrid = soln.grid.state(2,:);
% thGrid = soln.grid.state(3,:);
% uGrid = soln.grid.control;
% 
% % Plot the entire trajectory
% plot(x,y,'r-','LineWidth',3);
% 
% % Plot the grid points:
% plot(xGrid, yGrid, 'ko','MarkerSize',5,'LineWidth',3);
% 
% % Plot the start and end points:
% plot(x([1,end]), y([1,end]),'ks','MarkerSize',12,'LineWidth',3);
% 
% % Plot the state and control:
% figure(2); clf; 
% 
% subplot(2,2,1); hold on;
% plot(t,x);
% plot(tGrid,xGrid,'ko','MarkerSize',5,'LineWidth',3);
% ylabel('x');
% 
% subplot(2,2,3); hold on;
% plot(t,y);
% plot(tGrid,yGrid,'ko','MarkerSize',5,'LineWidth',3);
% ylabel('y');
% 
% subplot(2,2,2); hold on;
% plot(t,th);
% plot(tGrid,thGrid,'ko','MarkerSize',5,'LineWidth',3);
% ylabel('θ');
% 
% subplot(2,2,4); hold on;
% plot(tGrid,uGrid,'ko','MarkerSize',5,'LineWidth',3);
% plot(t,u);
% ylabel('u');