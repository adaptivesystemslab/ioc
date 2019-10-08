function model = createJumpModel_ioc_2D(JA, targ, jump, EKFCodePath)
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


% visu = 0; %for debugging/stepping through function
% 
% markerNames = fieldnames(trc.data);
% modelLinks = [];
% 
% %Scale to meters --> NO, HJC equation from referenced source is in mm
% 
% %First seconds just standing still, from this we can calculate the
% %world_to_base transform
% 
% crop1 = crop(1);
% crop2 = crop(2);
% 
% %Rotation calculate based on the mean of mid->front and right for y
% %Here we get the world to base transform to match Mocap Data
% front = mean((trc.data.HR(crop1:crop2,:)+trc.data.HL(crop1:crop2,:))/2);
% back = mean((trc.data.BR(crop1:crop2,:)+trc.data.BL(crop1:crop2,:))/2);
% mid = (front+back)/2;
% left = mean(trc.data.HL(crop1:crop2,:));
% right = mean(trc.data.HR(crop1:crop2,:));
% 
% %     P = [front(1:2) 0];
% %     Q = [mid(1:2) 0];
% %     R = [];
% 
% %This is rotation from world to body orieantation in mocap
% [~,R] = points2rot([front(1:2) 0],[mid(1:2) 0],[left(1:2) 0]);
% %R = eye(3);
% %     [x y z] = dcm2angle(R,'XYZ');
% 
% world_to_base = eye(4);
% world_to_base(1:3,1:3) = R'*roty(pi/2); %need the roty to keep joints correct
% world_to_base(1:3,4) = mid/1000;
% %R = eye(3);
% lasis = mean(R*trc.data.HL(crop1:crop2,:)',2);
% rasis = mean(R*trc.data.HR(crop1:crop2,:)',2);
% lback = mean(R*trc.data.BL(crop1:crop2,:)',2);
% rback = mean(R*trc.data.BR(crop1:crop2,:)',2);
% 
% %     figure(2); clf; hold on; grid on;
% %     plot3(lasis(1),lasis(2),lasis(3),'c*');
% %     plot3(rasis(1),rasis(2),rasis(3),'r*');
% %     plot3(lback(1),lback(2),lback(3),'g*');
% %     plot3(rback(1),rback(2),rback(3),'b*');
% %     xlabel('X');
% %     ylabel('Y');
% 
% %Distance between hips as calculated in 
% 
% %@ARTICLE{harrington2007prediction,
% %    author = {Harrington, ME and Zavatsky, AB and Lawson, SEM and Yuan, Z and Theologis,TN},
% %    title = {Prediction of the hip joint centre in adults, children, and patients
% %    with cerebral palsy based on magnetic resonance imaging},
% %    journal = {Journal of biomechanics},
% %    year = {2007},
% %    volume = {40}
% %}
% 
% %Pelvis Width calculated using x and y only
% PW = norm(lasis(1:2) - rasis(1:2));
% %Pelvis Depth calculated using x and y only
% PD = norm(abs((lasis(1:2)+rasis(1:2))/2 - (lback(1:2)+rback(1:2))/2));
% %Leg Length
% lankle = (mean(R*trc.data.ALLat(crop1:crop2,:)',2) + mean(R*trc.data.ALMed(crop1:crop2,:)',2))/2;
% rankle = (mean(R*trc.data.ARLat(crop1:crop2,:)',2) + mean(R*trc.data.ARMed(crop1:crop2,:)',2))/2;
% 
% LL = (norm(rasis-rankle) + norm(lasis-lankle))/2;
% 
% %So going from the middle the joint center is predicted in mm as
% x = -0.24*PD-9.9;
% y = -0.16*PW-0.04*LL-7.1;
% z = 0.28*PD+0.16*PW+7.9;
% % Where x is from middle to front, z is from middle to side, and y is
% % up
% 
% %Middle to RIGHT/LEFT hip JOINT CENTERS
% mid2RHC = ([x -z y]/1000)';
% mid2LHC = ([x z y]/1000)';
% 
% for m=3:numel(markerNames) %skip first 2 fields: frameNum, time
%        trc.data.(markerNames{m}) = (trc.data.(markerNames{m})')'./1000;
% end
% 
% 
% %Calculate joint centers
% lknee = (mean(trc.data.KLLat(crop1:crop2,:)) + mean(trc.data.KLMed(crop1:crop2,:)))/2;
% rknee = (mean(trc.data.KRLat(crop1:crop2,:)) + mean(trc.data.KRMed(crop1:crop2,:)))/2;
% lankle = (mean(trc.data.ALLat(crop1:crop2,:)) + mean(trc.data.ALMed(crop1:crop2,:)))/2;
% rankle = (mean(trc.data.ARLat(crop1:crop2,:)) + mean(trc.data.ARMed(crop1:crop2,:)))/2;
% ltoe = (mean(trc.data.FLLat(crop1:crop2,:)) + mean(trc.data.FLMed(crop1:crop2,:)))/2;
% rtoe = (mean(trc.data.FRLat(crop1:crop2,:)) + mean(trc.data.FRMed(crop1:crop2,:)))/2;
% shoulderCenterDrop = [0, 0, -S2sc];
% upperback = (mean(trc.data.SL(crop1:crop2,:)) + mean(trc.data.SR(crop1:crop2,:)))/2 + shoulderCenterDrop;
% lshoulder = mean(trc.data.SL(crop1:crop2,:)) + shoulderCenterDrop;
% rshoulder = mean(trc.data.SR(crop1:crop2,:)) + shoulderCenterDrop;
% lelbow = (mean(trc.data.ELLat(crop1:crop2,:)) + mean(trc.data.ELMed(crop1:crop2,:)))/2;
% relbow = (mean(trc.data.ERLat(crop1:crop2,:)) + mean(trc.data.ERMed(crop1:crop2,:)))/2;
% lwrist = (mean(trc.data.WLLat(crop1:crop2,:)) + mean(trc.data.WLMed(crop1:crop2,:)))/2;
% rwrist = (mean(trc.data.WRLat(crop1:crop2,:)) + mean(trc.data.WRMed(crop1:crop2,:)))/2;











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

% model = rlCModel('JumpModel_Eul.xml');
model = rlCModel('JumpModel_Eul_inertia_2D.xml');
transformNames = {model.transforms.name};
model.forwardPosition;

%     initPos = zeros(size(model.position));
initPos = zeros(numel(model.joints),1);

% because of 2D representation, just do L-pose for everyone
initPos(5) = pi/2; %shoulder elevation
initPos(6) = pi/2; %elbow flexion
initPos(7) = pi/2;
initPos(8) = pi/2; 


if(visu)
    vis = rlVisualizer('visModel',640,960);
    vis.addModel(model);
    vis.update();
end

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
modelLinks.spine = norm(upperback - FT(1:3,4)');
%     vis.update();
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'upperBack2lshoulder')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*(lshoulder - FT(1:3,4)')';
T(1:3,4) = FT(1:3,1:3)'*(modelLinks_ave.spine2shldr_L)';
model.transforms(trNum).t = T;
model.forwardPosition;
modelLinks.spine2shldr_L = (lshoulder - FT(1:3,4)'); % negative because extends down in global Z direction
%     vis.update();

trNum = find(ismember(transformNames,'upperBack2rshoulder')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*(rshoulder - FT(1:3,4)')';
T(1:3,4) = FT(1:3,1:3)'*(modelLinks_ave.spine2shldr_R)';
model.transforms(trNum).t = T;
model.forwardPosition;
modelLinks.spine2shldr_R = (rshoulder - FT(1:3,4)');
%     vis.update();


%% Pelvis
trNum = find(ismember(transformNames,'lAsis2Hip')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t; %mid_asis
% T(1:3,4) = mid2LHC;
T(1:3,4) = modelLinks_ave.base2hip_L;
model.transforms(trNum).t = T;
model.forwardPosition;
modelLinks.base2hip_L = mid2LHC';
%     vis.update();
% model = setInertialParam(model, trNum, gender, weight); % pelvis CoM only
% added with right hip link body

trNum = find(ismember(transformNames,'rAsis2Hip')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t; %mid_asis
% T(1:3,4) = mid2RHC;
T(1:3,4) = modelLinks_ave.base2hip_R;
model.transforms(trNum).t = T;
model.forwardPosition;
modelLinks.base2hip_R = mid2RHC';
%     vis.update();
model = setInertialParam(model, trNum, gender, weight);


%% Left Arm
trNum = find(ismember(transformNames,'lShoulder2Elbow')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(lshoulder - lelbow)]';
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.uparm_L]';
model.transforms(trNum).t = T;
model.forwardPosition;
modelLinks.uparm_L = -norm(lshoulder - lelbow);
%     vis.update();
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'lElbow2Wrist')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(lelbow - lwrist)]';
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.forearm_L]';
model.transforms(trNum).t = T;
model.forwardPosition;
modelLinks.forearm_L = -norm(lelbow - lwrist);
%     vis.update();
model = setInertialParam(model, trNum, gender, weight);


%% Right Arm
trNum = find(ismember(transformNames,'rShoulder2Elbow')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(rshoulder - relbow)]';
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.uparm_R]';
model.transforms(trNum).t = T;
model.forwardPosition;
modelLinks.uparm_R = -norm(rshoulder - relbow);
%     vis.update();
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'rElbow2Wrist')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(relbow - rwrist)]';
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.forearm_R]';
model.transforms(trNum).t = T;
model.forwardPosition;
modelLinks.forearm_R = -norm(relbow - rwrist);
%     vis.update();
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
modelLinks.thigh_L = -norm(T(1:3,4));
%     vis.update();
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'lKnee2Ankle')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(lknee - lankle)]';
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.shin_L]';
model.transforms(trNum).t = T;
model.forwardPosition;
modelLinks.shin_L = -norm(lknee - lankle);
%     vis.update();
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'lAnkle2Toe')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*(ltoe - FT(1:3,4)')';
T(1:3,4) = FT(1:3,1:3)'*(modelLinks_ave.foot_L)';
model.transforms(trNum).t = T;
model.forwardPosition;
modelLinks.foot_L = (ltoe - FT(1:3,4)');
%     vis.update();
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
modelLinks.thigh_R = -norm(T(1:3,4));
%     vis.update();
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'rKnee2Ankle')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(rknee - rankle)]';
T(1:3,4) = FT(1:3,1:3)'*[0,0,modelLinks_ave.shin_R]';
model.transforms(trNum).t = T;
model.forwardPosition;
modelLinks.shin_R = -norm(rknee - rankle);
%     vis.update();
model = setInertialParam(model, trNum, gender, weight);

trNum = find(ismember(transformNames,'rAnkle2Toe')==1);
T = model.transforms(trNum).t;
FT = model.transforms(trNum).frame_in.t;
% T(1:3,4) = FT(1:3,1:3)'*(rtoe - FT(1:3,4)')';
T(1:3,4) = FT(1:3,1:3)'*(modelLinks_ave.foot_R)';
model.transforms(trNum).t = T;
model.forwardPosition;
modelLinks.foot_R = (rtoe - FT(1:3,4)');
%     vis.update();
model = setInertialParam(model, trNum, gender, weight);

