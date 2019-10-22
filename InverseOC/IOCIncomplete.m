function IOCIncomplete(trialInfo,savePath)

%% -----------------load the model and other information----------------
% get the dynamics model of human body.
model = getModel(trialInfo);
% get basic information for IOC
trialInfo.numWeights = length(trialInfo.candidateFeatures);
trialInfo.numDofs = model.getModelDof();
%-----------------load the original trajectory data------------------
% load and generate the trajectory for the whole horizon
traj = model.loadData(trialInfo);
% Initialization of threshold and condition variables
dt = traj.trajT(2) - traj.trajT(1);
clc


%% -----------------cut the orginal trajectory manually----------------
% cutSavePath=fullfile('..\Data\IOC\CuttedData\Squat_IIT\',trialInfo.runName(1:6),'up');  
% cutSavePath=fullfile('..\Data\IOC\CuttedData\KneeHipFlexion\',trialInfo.runName(1:6),'up');  
cutSavePath=fullfile('..\Data\IOC\CuttedData\HipFlexion\',trialInfo.runName(1:6),'up');  
% cutSavePath=fullfile('..\Data\IOC\CuttedData\Squat\',trialInfo.runName(1:6),'up');  
% cutSavePath=fullfile('..\Data\IOC\CuttedData\Sit\',trialInfo.runName(1:6),'up');  

% ----------this section can also be commented if not needed----------
checkMkdir(cutSavePath) %check the folder exists
clearFlag=true; %clear or not the existing files
VisualizerCutter(model,traj,cutSavePath,clearFlag)
% VisualizerCutter(model,traj)


%% -----------------read the cropped trajectory------------------------
% read the cropped trajectory.
files=dir(fullfile(cutSavePath,'\*.mat'));
if(isempty(files)) return; end %#ok<SEPEX>
segIntervalIndex=load(fullfile(files(1).folder,files(1).name));
segIntervalIndex=segIntervalIndex.segIntervalIndex;
if(isempty(segIntervalIndex)) return; end %#ok<SEPEX>

%% ----------------perform ioc on each cutted traj---------------------
% Create IOC instance
ioc = IOCInstance(model, dt);
ioc.init(trialInfo);

% processing while increasing the window length
for k=1:size(segIntervalIndex,1)
    % take out the current segment
    segFrame=segIntervalIndex(k,1):segIntervalIndex(k,2);
    segTraj.trajX=traj.trajX(segFrame,:);
    segTraj.trajU=traj.trajU(segFrame,:);
    segTraj.trajT=traj.trajT(segFrame);
    segTraj.q=traj.q(segFrame,:);
    segTraj.dq=traj.dq(segFrame,:);
    segTraj.tau=traj.tau(segFrame,:);
    segTraj.frames=segFrame;
    
    
    
    % Initialize the recovery matrix [H1, H2];
    H1=[];
    H2=[];
    % IOC using the beginning segment of the trajectory
    for i=1:length(segFrame)
        % print
        clc
        sprintf('Data: %s', trialInfo.runName)
        sprintf('Prcentage: %.1f%%', i/length(segFrame)*100)
        currX=segTraj.trajX(i,:);
        currU=segTraj.trajU(i,:);
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
        Hhat = H;
        % solve for the weights from the current recvoery matirx
        [weights, ~] = computeWeights(Hhat, trialInfo.numWeights);
        if (isempty(trialInfo.gamma))
            trialInfo.gamma=4;
        end
        [rankInd, ~, ~] = validateCompletion(Hhat,trialInfo.gamma, trialInfo.delta, trialInfo.dimWeights);
        % storage for each step
        weightsVec(i,:)=weights;
        rankIndVec(i,:)=rankInd;
    end
    
    %store the results
    iocResults(k).weightsVec=weightsVec;
    iocResults(k).rankIndVec=rankIndVec;
    iocResults(k).traj=segTraj;
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
    segIntervalMat=[];
end
for i = 1:length(traj.trajT)
    %visualize
    mdl.position = traj.q(i, :);
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration();
    mdl.inverseDynamics();
    vis.update();
    clc
     model.model.joints.name
    rad2deg(traj.q(i,:))
    if (i==1)
        pause
    end
     pause;
    
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
            segIntervalMat(end+1,:)=interval
            indStart=[];
            indEnd=[];
        end
    end
    
    
    pause(0.01)
end


%save the segment
segIntervalIndex=[];
if cut && ~isempty(segIntervalMat)
    checkMkdir(cutSavePath) 
    
    if clearFlag
       delete(fullfile(cutSavePath,'*.mat')) %delete all the mat files under the folder
    end
    for k=1:size(segIntervalMat,1)
        %save
        t=segIntervalMat(k,1):segIntervalMat(k,2);
        
        if t(end)<length(traj.trajT)
            
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
                segIntervalIndex(end+1,:)=segIntervalMat(k,:);
%                 segTraj.q=traj.q(t,:);
%                 segTraj.dq=traj.dq(t,:);
%                 segTraj.tau=traj.tau(t,:);
%                 segTraj.trajT=traj.trajT(t,:);
%                 segTraj.trajU=traj.trajU(t,:);
%                 segTraj.trajX=traj.trajX(t,:);
%                 segTraj.frameInds=traj.frameInds(t)-traj.frameInds(t(1))+1;
% %                 segTraj.qLabels=traj.qLabels;
                savefile=strcat(cutSavePath,'/segIntervalIndex.mat');
                save(savefile,'segIntervalIndex')
                segIntervalIndex
                disp('Saved!')
            else
                continue;
            end
        end
        
    end
end

end




