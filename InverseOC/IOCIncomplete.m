function IOCIncomplete(trialInfo,savePath)

%-----------------load the model and other information----------------
% get the dynamics model of human body.
model = getModel(trialInfo);
% get basic information for IOC
trialInfo.numWeights = length(trialInfo.candidateFeatures);
trialInfo.numDofs = model.getModelDof();
%-----------------load the original trajectory data------------------
% load and generate the trajectory for the whole horizon
[q, dq, ddq, tau, trajT, trajU, trajX, traj] = model.loadData(trialInfo);
% Initialization of threshold and condition variables
dt = traj.trajT(2) - traj.trajT(1);

%-----------------cut the orginal trajectory manually----------------
% % set the savefule for the trimed trajectory
% cutSavePath=fullfile('..\Data\IOC\CuttedData\SquatMotion\',trialInfo.runName(1:11));   % for squatting motion
% cutSavePath=fullfile('..\Data\IOC\CuttedData\Jumping',trialInfo.runName(1:26));   % for the jumping motion
% %----------this section can also be commented if not needed----------
% checkMkdir(cutSavePath) %check the folder exists
% clearFlag=true; %clear or not the existing files
% VisualizerCutter(model,traj,cutSavePath,clearFlag)
VisualizerCutter(model,traj)
return 

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




