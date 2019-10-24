clc;clear all;close all;

% Add paths to directories with model definition and util functions
setPaths();

% Load ioc with information about each trial
configFile = jsondecode(fileread('../Data_json/PamelaConfig/DOC_Template2.json'));


% Create setup and define continuous functions for GPOPS2
currentTrial = configFile.Trials;
weights=[0, 1, 0]';


% Create dynamic model and IOC instance associated to it. This
% instance will be used to compute target features
dynModel = getModel(currentTrial);
dt=0.01;
ioc = IOCInstance(dynModel, dt);
ioc.init(currentTrial);

% Load the actual trajectory
traj=dynModel.loadData(configFile.Trials);
% read the cropped trajectory.
cutSavePath=fullfile('..\Data\IOC_Seg\CuttedData\KneeHipFlexion\Subj01\up'); 
files=dir(fullfile(cutSavePath,'\*.mat'));
if(isempty(files)) return; end %#ok<SEPEX>
segIntervalIndex=load(fullfile(files(1).folder,files(1).name));
segIntervalIndex=segIntervalIndex.segIntervalIndex;
if(isempty(segIntervalIndex)) return; end %#ok<SEPEX>
segFrame=segIntervalIndex(1,1):segIntervalIndex(1,2);
segTraj.trajX=traj.trajX(segFrame,:);
segTraj.trajU=traj.trajU(segFrame,:);
segTraj.trajT=traj.trajT(segFrame);
segTraj.q=traj.q(segFrame,:);
segTraj.dq=traj.dq(segFrame,:);
segTraj.tau=traj.tau(segFrame,:);
segTraj.frames=segFrame;



% Define optimization setup based on json entry
[q, dq, tau, time] = getSplineGuess(dynModel, currentTrial);

currentTrial.guess = struct('jointAngles', q, 'angularVelocities', dq,...
               'control', tau, 'time', time);
continuousFunc = @(input)humanDynContinuousFunc(input, weights, ioc);
setup = defineGPOPSOptimization(currentTrial, continuousFunc);

%  Solve problem using GPOPS2
output = gpops2(setup);

% Extract Solution and save it in data folder
solution = output.result.solution;
time = solution.phase.time;
state = solution.phase.state;
control = solution.phase.control;
% plot
figure(1)
clf
plot(time,state(:,1:2))
legend('q1','q2')
hold on
plot(segTraj.trajT-segTraj.trajT(1),segTraj.q,'.')
legend('actual_q1','actual_q2')
figure(2)
clf
plot(time,state(:,3:4))
figure(3)
clf
plot(time,control)
legend('tau1','tau2')








