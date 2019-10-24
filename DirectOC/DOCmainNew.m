clc;clear all;close all;

% Add paths to directories with model definition and util functions
setPaths();

% Load jin's ioc results
iocPath='../Data/JinIOCResults/IOC_KneeHipFlexion/Sub1/Subj01_KHEF_03CF.mat';
% iocPath='../Data/JinIOCResults/IOC_KneeHipFlexion/Sub1/Subj01_KHEF_04CF.mat';
% iocPath='../Data/JinIOCResults/IOC_KneeHipFlexion/Sub1/Subj01_KHEF_08CF.mat';
% iocPath='../Data/JinIOCResults/IOC_Squat_IIT/Sub1/Subj01_3DOF_3CF.mat';
% iocPath='../Data/JinIOCResults/IOC_Squat_IIT/Sub1/Subj01_3DOF_4CF.mat';
% iocPath='../Data/JinIOCResults/IOC_Squat_IIT/Sub1/Subj01_3DOF_8CF.mat';
load(iocPath)


% The original trajectory base directory
origDataBasePath='H:/data';
origDataSubPath=iocResultsSaver.origDataSubPath;
origDataPath=fullfile(origDataBasePath,origDataSubPath);

% Create setup the dynamics model and load the original data for doc
% results comparison
modelInfo.baseModel=iocResultsSaver.baseModel;
modelInfo.model=iocResultsSaver.model;
modelInfo.path=origDataPath;
modelInfo.runInds=[];
dynModel = getModel(modelInfo);
origTraj=dynModel.loadData(modelInfo);

% Create IOC instance to compute target features
dt=iocResultsSaver.dt;
ioc = IOCInstance(dynModel, dt);
featureInfo.candidateFeatures=iocResultsSaver.candidateFeatures;
ioc.init(featureInfo);

% Regenerate each segment of the trajectory by DOC, each loop is for one
% segment
for k=1:length(iocResultsSaver.iocResults)
    
    % DOC configuration setting
    iocResultsOfThisSeg=iocResultsSaver.iocResults(k);
    weights=iocResultsOfThisSeg.weights;
    config=struct;
    config.initialState.time=0;
    config.initialState.jointAngles=iocResultsOfThisSeg.initialJointAngle';
    config.initialState.angularVelocities=iocResultsOfThisSeg.initialJointVelocity';
    config.finalState.time=iocResultsOfThisSeg.finalTime-iocResultsOfThisSeg.initialTime; % watch here
    config.finalState.jointAngles=iocResultsOfThisSeg.finalJointAngle';
    config.finalState.angularVelocities=iocResultsOfThisSeg.finalJointVelocity';
    config.dt=dt;
    config.name=iocResultsSaver.iocConfig;
    config.bounds.minState=[-3, -3, -20, -20]';
    config.bounds.maxState=[3, 3, 20, 20]';
    config.bounds.minControl=[-1000, -1000]';
    config.bounds.maxControl=[1000, 1000]';
    
    % Define optimization setup based on json entry
    [q, dq, tau, time] = getSplineGuess(dynModel, config);
    config.guess = struct('jointAngles', q, 'angularVelocities', dq,...
        'control', tau, 'time', time);
    
    
    % Define dynamics function
    continuousFunc = @(input)humanDynContinuousFunc(input, weights, ioc);
    
    
    % Setup GPOPS-II configration
    setup = defineGPOPSOptimization(config, continuousFunc);
    %  Solve problem using GPOPS2
    output = gpops2(setup);
    
    % Extract Solution and save it in data folder
    solution = output.result.solution;
    time = solution.phase.time;
    state = solution.phase.state;
    control = solution.phase.control;
    
    
    % Plot the trajectory
    origSegJointAngles=origTraj.q(iocResultsOfThisSeg.trajSegIndex,:);
    origSegJointVelocity=origTraj.dq(iocResultsOfThisSeg.trajSegIndex,:);
    origSegTorque=origTraj.tau(iocResultsOfThisSeg.trajSegIndex,:);
    origSegTime=origTraj.trajT(iocResultsOfThisSeg.trajSegIndex)-origTraj.trajT(iocResultsOfThisSeg.trajSegIndex(1));
    
    figure(1)
    clf
    plot(time,state(:,1:2))
    legend('q1','q2')
    hold on
    plot(origSegTime,origSegJointAngles,'.')
    legend('orig_q1','orig_q2')
    
    figure(2)
    clf
    plot(time,state(:,3:4))
    legend('dq1','dq2')
    hold on
    plot(origSegTime,origSegJointVelocity,'.')
    legend('orig_dq1','orig_dq2')
    
    figure(2)
    clf
    plot(time,control(:,1:2))
    legend('tau1','tau2')
    hold on
    plot(origSegTime,origSegTorque,'.')
    legend('orig_tau1','orig_tau2')
    
end










