clear
clc
close all


setPaths();
iocPath='../Data/JinIOCResults/IOC_KneeHipFlexion/Sub1/Subj01_KHEF_03CF.mat';
load(iocPath)

% The original trajectory base directory
origDataBasePath='H:/data';
origDataSubPath=iocResultsSaver.origDataSubPath;
origDataPath=fullfile(origDataBasePath,origDataSubPath);

% Create setup the dynamics model and load the original data for doc
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

% Test new functions
iocResults=iocResultsSaver.iocResults(1);
weights=iocResults.weights;
initialState=[iocResults.initialJointAngle, iocResults.initialJointVelocity];

% Generate a initial control trajectory
config.initialState.time=0;
config.initialState.jointAngles=iocResults.initialJointAngle';
config.initialState.angularVelocities=iocResults.initialJointVelocity';
config.finalState.time=iocResults.finalTime-iocResults.initialTime; % watch here
config.finalState.jointAngles=iocResults.finalJointAngle';
config.finalState.angularVelocities=iocResults.finalJointVelocity';
config.dt=dt;
[qGuess, dqGuess, tauGuess, time] = getInitialGuess(dynModel, config,2);



% Compute the control grident
currControlTraj=tauGuess;
stop_thres=1e-1;
prevCost=inf;
count=1;
% Evaluate the current objective function
currCost=costEval(currControlTraj,initialState,ioc,weights);
fprintf("Initial cost is: %.4f\n", currCost);
while (abs(prevCost-currCost)>stop_thres)
    
    prevCost=currCost;
    % Compute the gradient on current control trajectory
    gradControl=trajGradient(currControlTraj,initialState,ioc,weights);
    
    % Update current control
    alpha=1/max(norm(gradControl,'fro'),100);
    nextControl=currControlTraj-alpha*gradControl;
    
    % Next loop
    currControlTraj=nextControl;
    count=count+1;
    [currCost,currStateTraj]=costEval(currControlTraj,initialState,ioc,weights);
    
    fprintf("Iteration %d, the current cost is: %.4f \n",count, currCost);
    
end


% Plot the trajectory
origSegJointAngles=origTraj.q(iocResults.trajSegIndex,:);
origSegJointVelocity=origTraj.dq(iocResults.trajSegIndex,:);
origSegTorque=origTraj.tau(iocResults.trajSegIndex,:);
origSegTime=origTraj.trajT(iocResults.trajSegIndex)-origTraj.trajT(iocResults.trajSegIndex(1));

figure(1)
clf
subplot(2,1,1);
plot(time,currStateTraj(:,1:size(currStateTraj,2)/2),'.-');
hold on
plot(origSegTime,origSegJointAngles,'.')
legend('q1','q2','origQ1','origQ2')




subplot(2,1,2);
plot(time,currControlTraj);
hold on
plot(origSegTime,origSegTorque,'.')
legend('tau1','tau2','origTau1','origTau2')


