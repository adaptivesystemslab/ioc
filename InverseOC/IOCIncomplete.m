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
% cutSavePath=fullfile('..\Data\IOC\CuttedData\SquatMotion\',trialInfo.runName(1:11));   % for squatting motion
cutSavePath=fullfile('..\Data\IOC\CuttedData\Jumping',trialInfo.runName(1:26));   % for the jumping motion
% %----------this section can also be commented if not needed----------
% checkMkdir(cutSavePath) %check the folder exists
% clearFlag=true; %clear or not the existing files
% VisualizerCutter(model,traj,cutSavePath,clearFlag)
% VisualizerCutter(model,traj)
% return 

%-----------------read the cropped trajectory------------------------
% read the cropped trajectory.
% if the squating motion
% files=dir(fullfile(cutSavePath,'\down\*.mat'));
% files=dir(fullfile(cutSavePath,'\down\*.mat'));

% if the jumping motion
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
        sprintf('Data: %s', trialInfo.runName)
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
%         Hhat = H;
        Hhat = H/norm(H,'fro');
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

switch trialInfo.baseModel
    case {'IIT','Jumping'}
        switch trialInfo.model
            case 'Jumping2D'
                load(trialInfo.path);   
                targNum = trialInfo.targNum;
                jumpNum = trialInfo.jumpNum;
                % Lin's code
%                 param.jump.takeoffFrame = JA.TOFrame(jumpNum,targNum);
%                 param.jump.landFrame = JA.LandFrame(jumpNum,targNum);
%                 param.jump.locationLand = JA.locationLand(12*(targNum-1) + jumpNum);
%                 param.jump.grade = JA.jumpGrades(12*(targNum-1) + jumpNum);
%                 param.jump.modelLinks = JA.modelLinks;
%                 param.jump.world2base = squeeze(JA.world2base((12*(targNum-1) + jumpNum),:,:));
%                 param.jump.bad2d = JA.bad2D(jumpNum,targNum);
%                 
%                 % Crop out initial and final calibration motions
%                 %                     takeoffFrames = 200; % 1 second before takeoff frame ...
%                 %                     landFrames = 300; % ... to 1.5 seconds after takeoff frame (~ 1 second after landing)
%                 %                     framesToUse = (param.jump.takeoffFrame-takeoffFrames):(param.jump.takeoffFrame+landFrames);
%                 
%                 takeoffFrames = 0; % 1 second before takeoff frame ...
%                 landFrames = 0; % ... to 1.5 seconds after takeoff frame (~ 1 second after landing)
%                 framesToUse = (param.jump.takeoffFrame-takeoffFrames):(param.jump.landFrame+landFrames);
%                 
%                 if(numel(framesToUse) < size(JA.targ(targNum).jump(jumpNum).data,1))
%                     fullDataAngles = JA.targ(targNum).jump(jumpNum).data(framesToUse,:);
%                 else % jump recording stops sooner than "landFrames" after TOFrame
%                     fullDataAngles = JA.targ(targNum).jump(jumpNum).data( framesToUse(1):end ,:);
%                     fullDataAngles = [fullDataAngles; repmat(fullDataAngles(end,:),(numel(framesToUse) - size(fullDataAngles,1)),1)]; % repeat last joint angle measurement for remainder of frames
%                 end
%                
                % Wanxin's modification
                fullDataAngles = JA.targ(targNum).jump(jumpNum).data;
                
                % keep only a subset of the joint angles
                qInds = [];
                allJointStr = {model.model.joints.name}';
                for indQ = 1:length(allJointStr)
                    qInds(indQ) = find(ismember(model.modelJointNameRemap, allJointStr{indQ}));
                end
                
                % also, negate the following joints since they're past
                % the flip
                qFlip = fullDataAngles;
                %                     jointsToFlip = {'rankle_jDorsiflexion', 'rknee_jExtension', 'rhip_jFlexion'};
                % %                     jointsToFlip = {'rankle_jDorsiflexion', 'rknee_jExtension', 'rhip_jFlexion', 'back_jFB', 'rjoint1'};
                %                     for indQ = 1:length(jointsToFlip)
                %                         qIndsFlip(indQ) = find(ismember(model.modelJointNameRemap, jointsToFlip{indQ}));
                %                         qFlip(:, qIndsFlip(indQ)) = -fullDataAngles(:, qIndsFlip(indQ));
                %                     end
                
                dt = 0.005;
                time = dt*(0:(size(qFlip, 1)-1));
                qRaw = qFlip(:, qInds);
                q = filter_dualpassBW(qRaw, 0.04, 0, 5);
                
                dqRaw = calcDerivVert(q, dt);
                dq = filter_dualpassBW(dqRaw, 0.04, 0, 5);
                %             dq = dqRaw;
                
                % don't filter ddq and tau to keep
                ddqRaw = calcDerivVert(dq, dt);
                %             ddq = filter_dualpassBW(ddqRaw, 0.04, 0, 5);
                ddq = ddqRaw;
                
                tauRaw = zeros(size(q));
                for indTime = 1:length(time) % recalc torque given redistributed masses
                    %                 model.updateState(q(indTime, :), dq(indTime, :));
                    tauRaw(indTime, :) = model.inverseDynamicsQDqDdq(q(indTime, :), dq(indTime, :), ddq(indTime, :));
                end
                
                %             tau = filter_dualpassBW(tauRaw, 0.04, 0, 5);
                tau = tauRaw;
                
                %             states = [q dq];
                states = encodeState(q, dq);
                control = tau;
                
                trajT = time';
                trajU = control;
                trajX = states;
                
                %output
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
                traj.frameInds=1:length(time);
                
                
            otherwise %case of "IIT"
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
%                 segTraj.qLabels=traj.qLabels;
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




