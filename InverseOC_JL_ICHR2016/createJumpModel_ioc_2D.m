function [model, initPos] = createJumpModel_ioc_2D(JA, targ, jump, EKFCodePath, modelCopyDOC)
%Calculate link lengths for Jumping Data model.
% uses frames between crop(1) and crop(2) of steady standing data


% TRC marker legend (center and right side):
% Neck          = back and base of neck (spine C7)
% SR            = shoulder right
% ERLat/ERMed   = elbow right lateral/medial
% WRLat/WRMed   = wrist right lateral/medial
% HR or HR  = right hip (anterior superior iliac spine)
% KRLat/KRMed   = knee right lateral/medial
% ARLat/ARMed   = ankle right lateral/medial
% FRLat         = foot right lateral (base of little toe)
% FRMed         = foot right medial (base of big toe)
% FRHeel        = foot right heel (back of calcaneous)


if(modelCopyDOC)
    model = rlCModel('JumpModel_Eul_inertia_2D_DOC.xml');
else
    model = rlCModel('JumpModel_Eul_inertia_2D.xml');
end


%% Create model file and attach markers
%     addpath(cd(cd('..')));
addpath(EKFCodePath);


%% Load model link lengths for participant
% mydir = pwd;
% idcs = strfind(mydir,'\');
% newdir = mydir(1:idcs(end)); 
% loadFilePath = [newdir 'results\Model_Info_JA_P' partNum '.mat'];
% load(loadFilePath);

modelLinks_ave = JA.modelLinks;

% [age,gender,height,weight,S2sc,calibPose] = getPartData(partNum);

gender = JA.partData.gender;
% height = JA.partData.height;
weight = JA.partData.weight;
% S2sc = JA.partData.S2sc;

% model = rlCModel('JumpModel_Eul_inertia_2D.xml');
transformNames = {model.transforms.name};
model.forwardPosition;

%     initPos = zeros(size(model.position));
initPos = zeros(numel(model.joints),1);

% because of 2D representation, just do L-pose for everyone
initPos(5) = pi/2; %shoulder elevation
initPos(6) = pi/2; %elbow flexion
initPos(7) = pi/2;
initPos(8) = pi/2; 


%World to pelvis center
trNum = find(ismember(transformNames,'world_to_base')==1);
model.transforms(trNum).t = squeeze(JA.world2base(12*(targ-1) + jump,:,:));
model.forwardPosition;

% vis.update();


frameList = {'back2Upperback','rAsis2Hip','lShoulder2Elbow','lElbow2Wrist',...
    'rShoulder2Elbow','rElbow2Wrist','lHip2Knee','lKnee2Ankle','lAnkle2Toe',...
    'rHip2Knee','rKnee2Ankle','rAnkle2Toe'};
for i = 1:numel(frameList)
    trNum = find(ismember(transformNames,frameList{i})==1);
    model = setInertialParam(model, trNum, gender, weight);
end

pause(0.01);



%% Neck and Shoulders
trNum = find(ismember(transformNames,'back2Upperback')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*([0,0,norm(upperback - FT(1:3,4)')])';
T(1:3,4) = FT(1:3,1:3)'*([0,0,modelLinks_ave.spine])';
model.transforms(trNum).t = T;
model.forwardPosition;
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'upperBack2lshoulder')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*(lshoulder - FT(1:3,4)')';
T(1:3,4) = FT(1:3,1:3)'*(modelLinks_ave.spine2shldr_L)';
model.transforms(trNum).t = T;
model.forwardPosition;

trNum = find(ismember(transformNames,'upperBack2rshoulder')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*(rshoulder - FT(1:3,4)')';
T(1:3,4) = FT(1:3,1:3)'*(modelLinks_ave.spine2shldr_R)';
model.transforms(trNum).t = T;
model.forwardPosition;


%% Pelvis
trNum = find(ismember(transformNames,'lAsis2Hip')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t; %mid_asis
% T(1:3,4) = mid2LHC;
T(1:3,4) = modelLinks_ave.base2hip_L;
model.transforms(trNum).t = T;
model.forwardPosition;
% model = setInertialParam(model, trNum, gender, weight); % pelvis CoM only
% added with right hip link body

trNum = find(ismember(transformNames,'rAsis2Hip')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t; %mid_asis
% T(1:3,4) = mid2RHC;
T(1:3,4) = modelLinks_ave.base2hip_R;
model.transforms(trNum).t = T;
model.forwardPosition;
model = setInertialParam(model, trNum, gender, weight);


%% Left Arm
trNum = find(ismember(transformNames,'lShoulder2Elbow')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(lshoulder - lelbow)]';
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.uparm_L]';
model.transforms(trNum).t = T;
model.forwardPosition;
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'lElbow2Wrist')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(lelbow - lwrist)]';
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.forearm_L]';
model.transforms(trNum).t = T;
model.forwardPosition;
model = setInertialParam(model, trNum, gender, weight);


%% Right Arm
trNum = find(ismember(transformNames,'rShoulder2Elbow')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(rshoulder - relbow)]';
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.uparm_R]';
model.transforms(trNum).t = T;
model.forwardPosition;
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'rElbow2Wrist')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(relbow - rwrist)]';
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.forearm_R]';
model.transforms(trNum).t = T;
model.forwardPosition;
model = setInertialParam(model, trNum, gender, weight);


%% Left Leg
trNum = find(ismember(transformNames,'lHip2Knee')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*(lknee - FT(1:3,4)')';
% % take norm of above, offset knee straight down, directly below hip
% T(1:3,4) = [0,0,-norm(T(1:3,4))];
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.thigh_L]';
model.transforms(trNum).t = T;
model.forwardPosition;
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'lKnee2Ankle')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(lknee - lankle)]';
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.shin_L]';
model.transforms(trNum).t = T;
model.forwardPosition;
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'lAnkle2Toe')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*(ltoe - FT(1:3,4)')';
T(1:3,4) = FT(1:3,1:3)'*(modelLinks_ave.foot_L)';
model.transforms(trNum).t = T;
model.forwardPosition;
model = setInertialParam(model, trNum, gender, weight);


%% Right Leg
trNum = find(ismember(transformNames,'rHip2Knee')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*(rknee - FT(1:3,4)')';
% % take norm of above, offset knee straight down, directly below hip
% T(1:3,4) = [0,0,norm(T(1:3,4))];
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.thigh_R]';
model.transforms(trNum).t = T;
model.forwardPosition;
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'rKnee2Ankle')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(rknee - rankle)]';
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.shin_R]';
model.transforms(trNum).t = T;
model.forwardPosition;
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'rAnkle2Toe')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*(rtoe - FT(1:3,4)')';
T(1:3,4) = FT(1:3,1:3)'*(modelLinks_ave.foot_R)';
model.transforms(trNum).t = T;
model.forwardPosition;
model = setInertialParam(model, trNum, gender, weight);

