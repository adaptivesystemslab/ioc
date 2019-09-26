function IOCIncomplete(trialInfo,savePath)

%-----------------load the model and other information----------------
% get the dynamics model of human body.
model = getModel(trialInfo);
% get basic information for IOC
trialInfo.numWeights = length(trialInfo.candidateFeatures);
trialInfo.numDofs = model.totalActuatedJoints;
%-----------------load the original trajectory data------------------
% load and generate the trajectory for the whole horizon
traj= loadData(trialInfo, model);
% Initialization of threshold and condition variables
dt = traj.trajT(2) - traj.trajT(1);

%-----------------cut the orginal trajectory manually----------------
% % set the savefule for the trimed trajectory
cutSavePath=fullfile('..\Data\IOC\CuttedData\',trialInfo.runName(1:11));
% %----------this section can also be commented if not needed----------
% checkMkdir(cutSavePath) %check the folder exists
% clearFlag=true; %clear or not the existing files
% VisualizerCutter(model,traj,cutSavePath,clearFlag)
% return 

%-----------------read the cropped trajectory------------------------
% read the cropped trajectory.
% files=dir(fullfile(cutSavePath,'\down\*.mat'));
files=dir(fullfile(cutSavePath,'\*.mat'));
for k=1:length(files)
    dataPath=fullfile(files(k).folder,files(k).name);
    trajSet(k)=load(dataPath);
end


%----------------perform ioc on each cutted traj---------------------
% Create IOC instance
ioc = IOCInstanceNew(model, dt);
ioc.init(trialInfo);
% initialize the result vector

for k=1:length(trajSet)
    traj=trajSet(k).segTraj;   
    % visualize the trajectory first
%     VisualizerCutter(model,traj)
    
    
    % Initialize the recovery matrix [H1, H2];
    H1=[];
    H2=[];
    % IOC using the beginning segment of the trajectory
    rankIndVec=[];
    weightsVec=[];
    for i=1:traj.frameInds(end)
        % print
        clc
        sprintf('Prcentage: %.1f%%', i/traj.frameInds(end)*100)
        currX=traj.trajX(i,:);
        currU=traj.trajU(i,:);
        %compute the deriviatives(jacobian matrix) for the system
        Jx=ioc.getDynamicsStateDerivatives(currX,currU);
        Ju=ioc.getDynamicsControlDerivatives(currX,currU);
        %compute the deriviatives(jacobina matrix) for the features
        Px=ioc.getFeaturesStateDerivatives(currX,currU);
        Pu=ioc.getFeaturesControlDerivatives(currX,currU);
        %take transpose for all derivatives
        Jx=Jx';
        Ju=Ju';
        Px=Px';
        Pu=Pu';
        % assemble the recovery matrix
        [H1,H2]=assembleH1H2(Jx, Ju, Px, Pu, H1,H2);
        H=[H1, H2];
        % check the status of the recovery matrix
        Hhat = H;
        % solve for the weights from the current recvoery matirx
        [weights, ~] = computeWeights(Hhat, trialInfo.numWeights);
        % solve for the rank index from the current recovery matrix
        if (isempty(trialInfo.gamma))
            trialInfo.gamma=4;
        end
        [rankInd, completed, errorCode] = validateCompletion(Hhat,trialInfo.gamma, trialInfo.delta, trialInfo.dimWeights);
        % storage for each step
        weightsVec(end+1,:)=weights;
        rankIndVec(end+1,:)=rankInd;
    end
    
    %store the results
    iocResults(k).weightsVec=weightsVec;
    iocResults(k).rankIndVec=rankIndVec;
    iocResults(k).traj=traj;
    for ii=1:length(ioc.features)
        iocResults(k).featureLabels{ii}=char(ioc.features(ii).feature);
    end
    
    % plot the results
    figure(1)
    clf
    for r=1:length(iocResults(end).featureLabels)
        subplot(length(iocResults(end).featureLabels),1,r)
        plot(iocResults(end).weightsVec(:,r));
        ylabel(iocResults(end).featureLabels{r})
        ylim([0 1])
    end
    drawnow
    
%     figure(2)
%     clf
%     plot(iocResults(end).rankIndVec);
%     drawnow

disp('pause')
pause;

end




%----------------------------save the results-----------------------------
% save('iocResults.mat','iocResults')


end

function traj = loadData(trialInfo, model)
% Load state, control and time trajectories to be analyzed.
load(trialInfo.path);

% keep only the joint angles corresponding
qInds = [];
allJointStr = {model.model.joints.name}';

for indQ = 1:length(allJointStr)
    qInds(indQ) = find(ismember(saveVar.jointLabels, allJointStr{indQ}));
end


if ~isempty(trialInfo.runInds)
    frameInds = trialInfo.runInds(1):trialInfo.runInds(2);
else
    frameInds = 1:length(saveVar.time);
end

if max(frameInds) > length(saveVar.time)
    frameInds = frameInds(1):length(saveVar.time);
end

time = saveVar.time;
qRaw = saveVar.jointAngle.array(:, qInds);
%this should be noted because the raw q data is in the order of  'joint_hip', 'joint_knee', 'joint_ankle'
%also the joint_ankle is compesated by an offset pi/2 to make the data into
%vertical.
% qRaw=flip(qRaw,2)+[pi/2*ones(size(qRaw,1),1), zeros(size(qRaw,1),1), zeros(size(qRaw,1),1) ];
% qRaw=flip(qRaw,2);
q = filter_dualpassBW(qRaw, 0.04, 0, 5);

dqRaw = calcDerivVert(q, saveVar.dt);
dq = filter_dualpassBW(dqRaw, 0.04, 0, 5);

% don't filter ddq and tau to keep
ddqRaw = calcDerivVert(dq, saveVar.dt);
ddq = ddqRaw;

tauRaw = zeros(size(q));
for indTime = 1:length(time) % recalc torque given redistributed masses
    tauRaw(indTime, :) = model.inverseDynamicsQDqDdq(q(indTime, :), dq(indTime, :), ddq(indTime, :));
end

tau = tauRaw;
states = encodeState(q, dq);
control = tau;
trajT = time;
trajU = control;
trajX = states;

% outputs
traj.q=q;
traj.dq=dq;
traj.ddq=ddq;
traj.tau=tau;
traj.states=states;
traj.control=control;
traj.time=time';
traj.trajT=trajT;
traj.trajU=trajU;
traj.trajX=trajX;
traj.frameInds=frameInds;
traj.qLabels={'joint_hip', 'joint_knee', 'joint_ankle'};



end

function Subvisualizer(traj)

if(isempty(traj))
    return
end
trajq=traj.q;

% define the link for the ankle, knee, and hip.
Lankle=1;
Lknee=1;
Lhip=1;

% drawing
figure(1)
for i=traj.frameInds
    q_hip=trajq(i,1);
    q_knee=trajq(i,2);
    q_ankle=trajq(i,3)+pi/2;
    Pankle=[0, 0];
    Pknee=Pankle+[Lankle*cos(q_ankle), Lankle*sin(q_ankle)];
    Phip=Pknee+[Lknee*cos(q_ankle+q_knee), Lknee*sin(q_ankle+q_knee)];
    Phead=Phip+[Lhip*cos(q_hip+q_knee+q_ankle), Lhip*sin(q_hip+q_knee+q_ankle)];
    %plot
    clf
    axis([-5, 5, 0, 5]);
    hold on
    plot([Pankle(1), Pknee(1)],[Pankle(2), Pknee(2)]);
    text(Pankle(1),Pankle(2),'ankle')
    plot([Pknee(1),Phip(1)],[Pknee(2),Phip(2)]);
    text(Pknee(1),Pknee(2),'knee')
    plot([Phip(1),Phead(1)],[Phip(2),Phead(2)]);
    text(Phip(1),Phip(2),'hip')
    text(Phead(1),Phead(2),'head')
    hold off 
    pause(0.001)
end

end

function VisualizerCutter(model,traj,cutSavePath,clearFlag)

% check if there is any need to cut
switch nargin
    case {3,4}
        cut=true;
    otherwise
        cut=false;
end


mdl = model.model;
vis = rlVisualizer('vis',960,640);
mdl.forwardPosition();
vis.addModel(mdl);
vis.addMarker('x-axis', [1 0 0], [0.2 0.2 0.2 1]);
vis.addMarker('y-axis', [0 1 0], [0.2 0.2 0.2 1]);
vis.addMarker('z-axis', [0 0 1], [0.2 0.2 0.2 1]);
vis.update();

% If cut needed, initialize the croping starting point and ending point,
% and delete all the filer under the savefolder
if cut
    indStart=[];
    indEnd=[];
    indMat=[];
end
for i = 1:length(traj.trajT)
    %visualizer
    mdl.position = traj.q(i, :);
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration();
    mdl.inverseDynamics();
    vis.update();
    
    %% if cut needed
    if cut
        %get the key response: rightarrow key means crop starting point;
        %leftarrow key means crop ending point. In between, the corresponding
        %data is saved
        rsp=getkey(0.01);
        if strcmpi(rsp.keyName,'RightArrow')
            indStart=i;
        end
        if strcmpi(rsp.keyName,'LeftArrow')
            indEnd=i;
        end
        if strcmpi(rsp.keyName,'s')
            break;
        end
        if ~isempty(indStart)&& ~isempty(indEnd) && indEnd>indStart
            interval=[indStart, indEnd];
            indMat(end+1,:)=interval
            indStart=[];
            indEnd=[];
        end
    end
    
    
    pause(0.01)
end


%save the segment
if cut && ~isempty(indMat)
    checkMkdir(cutSavePath) 
    
    if clearFlag
       delete(fullfile(cutSavePath,'*.mat')) %delete all the mat files under the folder
    end
    for k=1:size(indMat,1)
        %save
        t=indMat(k,1):indMat(k,2);
        
        if t(end)<traj.frameInds(end)
            
            % animate it again for confirmation
            clc
            prompt = 'Please confirm this segment. Do you want save (Y/N) or replay (R)?: ';
            str=[];
            while ~strcmpi(str,'y') && ~strcmpi(str,'n')
                for i=t
                    %visualizer
                    mdl.position = traj.q(i, :);
                    mdl.forwardPosition();
                    mdl.forwardVelocity();
                    mdl.forwardAcceleration();
                    mdl.inverseDynamics();
                    vis.update();
                    pause(0.01);
                end
                str = input(prompt,'s');
            end
            if (strcmpi(str,'y')&&t(2))
                segTraj.q=traj.q(t,:);
                segTraj.dq=traj.dq(t,:);
                segTraj.ddq=traj.ddq(t,:);
                segTraj.tau=traj.tau(t,:);
                segTraj.states=traj.states(t,:);
                segTraj.control=traj.control(t,:);
                segTraj.time=traj.time(t);
                segTraj.trajT=traj.trajT(t,:);
                segTraj.trajU=traj.trajU(t,:);
                segTraj.trajX=traj.trajX(t,:);
                segTraj.frameInds=traj.frameInds(t)-traj.frameInds(t(1))+1;
                segTraj.qLabels=traj.qLabels;
                savefile=strcat(cutSavePath,'\Segment_',num2str(traj.trajT(t(1))),'_',num2str(traj.trajT(t(end))),'.mat');
                save(savefile,'segTraj')
                disp('Saved!')
            else
                continue;
            end
        end
        
    end
end

end




