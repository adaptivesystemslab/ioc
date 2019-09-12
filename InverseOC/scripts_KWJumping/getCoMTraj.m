function CoMTraj = getCoMTraj(JA,targNum,jumpNum,EKFCodePath)
% Calculates CoM 3D position throughout entire jump trajectory
% EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';
addpath(EKFCodePath);

visualize = 0; % debugging

jointTraj = JA.targAlign(targNum).jump(jumpNum).data;
dataLength = size(jointTraj,1);
currJump = 12*(targNum-1) + jumpNum;

[~,~,height,weight,~,~] = getPartData(JA.partNum);

if(numel(JA.modelLinks.spine)==1)
    use_const_model_links = 1;
else
    use_const_model_links = 0;
end


% Segments (11 total): head/torso/pelvis, L/R up-arm, L/R forearm, L/R thigh, 
% L/R shin, L/R foot
segFrames = {'back2Upperback','lShoulder2Elbow','rShoulder2Elbow','lElbow2Wrist','rElbow2Wrist',...
    'lHip2Knee','rHip2Knee','lKnee2Ankle','rKnee2Ankle','lAnkle2Toe','rAnkle2Toe'};
segModelFrames = {'spine','uparm_L','uparm_R','forearm_L','forearm_R',...
    'thigh_L','thigh_R','shin_L','shin_R','foot_L','foot_R'};
segWeights = double(weight).*[0.578, 0.028, 0.028, 0.022, 0.022, ...
    0.1, 0.1, 0.0465, 0.0465, 0.0145*1.25, 0.0145*1.25]; % NOTE: foot segment weights "*1.25" to approximate weight of shoe

if(use_const_model_links)
    base2hipZ = abs(JA.modelLinks.base2hip_L(3) + JA.modelLinks.base2hip_R(3))/2;
    segProximalDist = [((0.47*(double(height)/100) + base2hipZ)*0.396 - base2hipZ)/JA.modelLinks.spine,... %CoM vertical distance from base frame (up along spine)
        0.436, 0.436, 0.686, 0.686, 0.433, 0.433, 0.433, 0.433, 0.429, 0.429]; %CoM proximal distances of segments 2 to 10 
else
    base2hipZ = abs(JA.modelLinks.base2hip_L(currJump,3) + JA.modelLinks.base2hip_R(currJump,3))/2;
    segProximalDist = [((0.47*(double(height)/100) + base2hipZ)*0.396 - base2hipZ)/JA.modelLinks.spine(currJump),... %CoM vertical distance from base frame (up along spine)
        0.436, 0.436, 0.686, 0.686, 0.433, 0.433, 0.433, 0.433, 0.429, 0.429]; %CoM proximal distances of segments 2 to 10 
end
% NOTE: forearm distance "0.686" = combined CoM location of forearm and 
% hand, expressed as % distance of forearm length only, from elbow

CoMTraj = zeros(dataLength,3);


%% Make model object, set link lengths
addpath(EKFCodePath);
model = rlCModel('JumpModel_Eul.xml');
transformNames = {model.transforms.name};
model.forwardPosition;

if(visualize)
    vis = rlVisualizer('vis', 640, 960);
    vis.addModel(model);
    vis.update;
end

modelLinkFrames = {'back2Upperback','upperBack2lshoulder','upperBack2rshoulder',...
    'lAsis2Hip','rAsis2Hip','lShoulder2Elbow','lElbow2Wrist','rShoulder2Elbow','rElbow2Wrist',...
    'lHip2Knee','lKnee2Ankle','lAnkle2Toe','rHip2Knee','rKnee2Ankle','rAnkle2Toe'};
% modelLinkNames = {'spine','spine2shldr_L','spine2shldr_R','base2hip_L','base2hip_R',...
%     'uparm_L','forearm_L','uparm_R','forearm_R',...
%     'thigh_L','shin_L','foot_L','thigh_R','shin_R','foot_R'};
modelLinkNames = fieldnames(JA.modelLinks);

trNum = find(ismember(transformNames,'world_to_base')==1);
model.transforms(trNum).t = JA.world2base(currJump,:,:);
model.forwardPosition;


for link = 1:numel(modelLinkFrames)
%     if(~contains(modelLinkNames{link},'world2base'))
        trNum = find(ismember(transformNames,modelLinkFrames{link})==1);
        T = model.transforms(trNum).t;
        FT = model.transforms(trNum).frame_in.t;
        
        if(use_const_model_links)
            if(numel(JA.modelLinks.(modelLinkNames{link}))==1) %just Z-global direction
                linkOffset = [0,0,JA.modelLinks.(modelLinkNames{link})];
            else
                linkOffset = JA.modelLinks.(modelLinkNames{link});
            end
        else
            if(numel(JA.modelLinks.(modelLinkNames{link})(currJump,:))==1) %just Z-global direction
                linkOffset = [0,0,JA.modelLinks.(modelLinkNames{link})(currJump)];
            else
                linkOffset = JA.modelLinks.(modelLinkNames{link})(currJump,:);
            end
        end
        
        T(1:3,4) = FT(1:3,1:3)'*(linkOffset');
        model.transforms(trNum).t = T;
        model.forwardPosition;

%         if(visualize)
%             vis.update;
%         end
%     end
end




%% Find full body CoM at each timestep of jump 
for i = 1:dataLength
    
    segCoM = zeros(3,numel(segFrames));
    model.position = jointTraj(i,:); % set model to current pos
    model.forwardPosition;
    
    % find pos of CoM of each body segment
    for seg = 1:numel(segFrames)
        trNum = find(ismember(transformNames,segFrames{seg})==1);
        T = model.transforms(trNum).t;
        FT = model.transforms(trNum).frame_in.t;
        inFrameCoM = segProximalDist(seg)*T(1:3,4);
        inFrameCoM_T = FT*[T(1:3,1:3),inFrameCoM;[0,0,0,1]];
        segCoM(:,seg) = inFrameCoM_T(1:3,4);
        
        if(visualize)
            vis.addMarker(segFrames{seg},segCoM(:,seg));
%             vis.update;
        end
        
    end
    
    % calculate full body CoM
    CoMTraj(i,:) = segCoM*(segWeights'/sum(segWeights));
    
    if(visualize)
        vis.addMarker('FB CoM',CoMTraj(i,:));
        vis.update;
    end
end






