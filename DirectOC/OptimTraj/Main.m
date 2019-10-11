%% First version of DOC code for RL model and Recovery matrix IOC

% clear all;
clc;
tic;

%% Set up internal parameters
% Add paths to directories with model definition and util functions
setPaths();

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                Read trial data and weights                              %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

status = 0;

% Create and/or look for folder where solutions are going to be saved
currentDate = erase(string(datetime("today")),"-");
savePath = sprintf("../../Data/OC/%s/", currentDate);

if ~exist(savePath, 'dir')
   status = mkdir(char(savePath)); 
end

frameInitPose = 0;
frameFinalPose = 100;

timeInit = 0;
timeFinal = 1;

if status
    % Load data and weights with hardcoded paths for now
    motion = load('/home/pcarreno/Documents/WaterlooPostDoc/Code/RecoveryMatrixIOC/outputData/Subj1_3DOF_3CF/data_000001_000999.mat');
    iocOutput = load('/home/pcarreno/Documents/WaterlooPostDoc/Code/RecoveryMatrixIOC/outputData/Subj1_3DOF_3CF/weights_000001_000999.mat');
    weightTrajectory = load('/home/pcarreno/Documents/WaterlooPostDoc/Code/RecoveryMatrixIOC/outputData/Subj1_3DOF_3CF/mat_results_cumulativeAllPass_result01_Subj1_3DOF_3CF_3DOF_3CF_Forward');
    
    % Get information about model and features
    trialInfo = iocOutput.trialInfo;
    
    % Create model and instantiate ioc instance
    model = getModel(trialInfo);
    dt = motion.dt;
    ioc = IOCInstance(model, dt);

    % Get bounds on state and control from motion file
    [state_bounds, control_bounds] = getBounds(motion);
    
    % Get initial and final state. Generate initial guess from spline
    initialState = [motion.q(frameInitPose,:), motion.dq(frameInitPose,:)];
    finalState = [motion.q(frameFinalPose,:), motion.dq(frameFinalPose,:)];
    
    guess = guessFromSpline(model, initialState, finalState);
    
    % Get weights to be used during objective computation. 
    % For now I assume weights correspond to avg recovery at finalState frame
    weights = weightTrajectory(frameFinalPose,:);
    
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                      Set up function handles                            %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,x,u)( Dynamics(x,u,ioc) );
problem.func.pathObj = @(t,x,u)( Objective(t,x,u,weights,ioc) );


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                 Set up bounds on state and control                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

numDof = length(currentTrial.initialState.jointAngles); 

problem.bounds.initialTime.low = currentTrial.initialState.time;
problem.bounds.initialTime.upp = currentTrial.initialState.time;
problem.bounds.finalTime.low = currentTrial.finalState.time;
problem.bounds.finalTime.upp = currentTrial.finalState.time;

problem.bounds.state.low = [arrayfun(@(x) deg2rad(x), currentTrial.bounds.minState(1:numDof))', currentTrial.bounds.minState(numDof+1:end)']';
problem.bounds.state.upp = [arrayfun(@(x) deg2rad(x), currentTrial.bounds.maxState(1:numDof))', currentTrial.bounds.maxState(numDof+1:end)']';

problem.bounds.initialState.low = [arrayfun(@(x) deg2rad(x), currentTrial.initialState.jointAngles)', currentTrial.initialState.angularVelocities']';
problem.bounds.initialState.upp = [arrayfun(@(x) deg2rad(x), currentTrial.initialState.jointAngles)', currentTrial.initialState.angularVelocities']';

problem.bounds.finalState.low = [arrayfun(@(x) deg2rad(x), currentTrial.finalState.jointAngles)', currentTrial.finalState.angularVelocities']';
problem.bounds.finalState.upp = [arrayfun(@(x) deg2rad(x), currentTrial.finalState.jointAngles)', currentTrial.finalState.angularVelocities']';

problem.bounds.control.low = currentTrial.bounds.minControl;
problem.bounds.control.upp = currentTrial.bounds.maxControl;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                 Initialize trajectory with guess                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

anglesGuess = []; velocityGuess = []; controlGuess = [];
for i=1:length(currentTrial.guess)
    temp = currentTrial.guess(i);
    anglesGuess = [anglesGuess, arrayfun(@(x) deg2rad(x), temp.jointAngles)];
    velocityGuess = [velocityGuess, temp.angularVelocities];
    controlGuess = [controlGuess, temp.control];    
end
    
if isfield(currentTrial.guess(1), 'time')
    timeGuess = currentTrial.guess(1).time';
else
    timeGuess = [currentTrial.initialState.time; currentTrial.finalState.time];
end
    
problem.guess.time = timeGuess';
problem.guess.state = [anglesGuess, velocityGuess]';
problem.guess.control = controlGuess';


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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Display the solution                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

figure(1); clf; hold on;

drawHills(xBnd,yBnd);

t = linspace(soln.grid.time(1), soln.grid.time(end), 150);
z = soln.interp.state(t);
x = z(1,:);
y = z(2,:);
th = z(3,:);
u = soln.interp.control(t);

tGrid = soln.grid.time;
xGrid = soln.grid.state(1,:);
yGrid = soln.grid.state(2,:);
thGrid = soln.grid.state(3,:);
uGrid = soln.grid.control;

% Plot the entire trajectory
plot(x,y,'r-','LineWidth',3);

% Plot the grid points:
plot(xGrid, yGrid, 'ko','MarkerSize',5,'LineWidth',3);

% Plot the start and end points:
plot(x([1,end]), y([1,end]),'ks','MarkerSize',12,'LineWidth',3);

% Plot the state and control:
figure(2); clf; 

subplot(2,2,1); hold on;
plot(t,x);
plot(tGrid,xGrid,'ko','MarkerSize',5,'LineWidth',3);
ylabel('x');

subplot(2,2,3); hold on;
plot(t,y);
plot(tGrid,yGrid,'ko','MarkerSize',5,'LineWidth',3);
ylabel('y');

subplot(2,2,2); hold on;
plot(t,th);
plot(tGrid,thGrid,'ko','MarkerSize',5,'LineWidth',3);
ylabel('θ');

subplot(2,2,4); hold on;
plot(tGrid,uGrid,'ko','MarkerSize',5,'LineWidth',3);
plot(t,u);
ylabel('u');
