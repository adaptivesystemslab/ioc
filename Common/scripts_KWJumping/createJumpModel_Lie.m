function [model, trc, initPos_lie] = createJumpModel(trc, attachMarkers, crop, S2sc, EKFCodePath)
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


    visu = 0; %for debugging/stepping through function
    
    model_transforms = [];
    
    markerNames = fieldnames(trc.data);

    %Scale to meters --> NO, HJC equation from referenced source is in mm
    
    %First seconds just standing still, from this we can calculate the
    %world_to_base transform
    
    crop1 = crop(1);
    crop2 = crop(2);
    
    %Rotation calculate based on the mean of mid->front and right for y
    %Here we get the world to base transform to match Mocap Data
    front = mean((trc.data.HR(crop1:crop2,:)+trc.data.HL(crop1:crop2,:))/2);
    back = mean((trc.data.BR(crop1:crop2,:)+trc.data.BL(crop1:crop2,:))/2);
    mid = (front+back)/2;
    left = mean(trc.data.HL(crop1:crop2,:));
    right = mean(trc.data.HR(crop1:crop2,:));
    
%     P = [front(1:2) 0];
%     Q = [mid(1:2) 0];
%     R = [];
    
    %This is rotation from world to body orieantation in mocap
    [~,R] = points2rot([front(1:2) 0],[mid(1:2) 0],[left(1:2) 0]);
    %R = eye(3);
%     [x y z] = dcm2angle(R,'XYZ');

    world_to_base = eye(4);
    world_to_base(1:3,1:3) = R'*roty(pi/2); %need the roty to keep joints correct
    world_to_base(1:3,4) = mid/1000;
    %R = eye(3);
    lasis = mean(R*trc.data.HL(crop1:crop2,:)',2);
    rasis = mean(R*trc.data.HR(crop1:crop2,:)',2);
    lback = mean(R*trc.data.BL(crop1:crop2,:)',2);
    rback = mean(R*trc.data.BR(crop1:crop2,:)',2);

%     figure(2); clf; hold on; grid on;
%     plot3(lasis(1),lasis(2),lasis(3),'c*');
%     plot3(rasis(1),rasis(2),rasis(3),'r*');
%     plot3(lback(1),lback(2),lback(3),'g*');
%     plot3(rback(1),rback(2),rback(3),'b*');
%     xlabel('X');
%     ylabel('Y');

    %Distance between hips as calculated in 
    
    %@ARTICLE{harrington2007prediction,
    %    author = {Harrington, ME and Zavatsky, AB and Lawson, SEM and Yuan, Z and Theologis,TN},
    %    title = {Prediction of the hip joint centre in adults, children, and patients
    %    with cerebral palsy based on magnetic resonance imaging},
    %    journal = {Journal of biomechanics},
    %    year = {2007},
    %    volume = {40}
    %}
    
    %Pelvis Width calculated using x and y only
    PW = norm(lasis(1:2) - rasis(1:2));
    %Pelvis Depth calculated using x and y only
    PD = norm(abs((lasis(1:2)+rasis(1:2))/2 - (lback(1:2)+rback(1:2))/2));
    %Leg Length
    lankle = (mean(R*trc.data.ALLat(crop1:crop2,:)',2) + mean(R*trc.data.ALMed(crop1:crop2,:)',2))/2;
    rankle = (mean(R*trc.data.ARLat(crop1:crop2,:)',2) + mean(R*trc.data.ARMed(crop1:crop2,:)',2))/2;
    
    LL = (norm(rasis-rankle) + norm(lasis-lankle))/2;
    
    %So going from the middle the joint center is predicted in mm as
    x = -0.24*PD-9.9;
    y = -0.16*PW-0.04*LL-7.1;
    z = 0.28*PD+0.16*PW+7.9;
    % Where x is from middle to front, z is from middle to side, and y is
    % up

    %Middle to RIGHT/LEFT hip JOINT CENTERS
    mid2RHC = ([x -z y]/1000)';
    mid2LHC = ([x z y]/1000)';
    
    for m=3:numel(markerNames) %skip first 2 fields: frameNum, time
           trc.data.(markerNames{m}) = (trc.data.(markerNames{m})')'./1000;
    end
    
    
    %Calculate joint centers
    lknee = (mean(trc.data.KLLat(crop1:crop2,:)) + mean(trc.data.KLMed(crop1:crop2,:)))/2;
    rknee = (mean(trc.data.KRLat(crop1:crop2,:)) + mean(trc.data.KRMed(crop1:crop2,:)))/2;
    lankle = (mean(trc.data.ALLat(crop1:crop2,:)) + mean(trc.data.ALMed(crop1:crop2,:)))/2;
    rankle = (mean(trc.data.ARLat(crop1:crop2,:)) + mean(trc.data.ARMed(crop1:crop2,:)))/2;
    ltoe = (mean(trc.data.FLLat(crop1:crop2,:)) + mean(trc.data.FLMed(crop1:crop2,:)))/2;
    rtoe = (mean(trc.data.FRLat(crop1:crop2,:)) + mean(trc.data.FRMed(crop1:crop2,:)))/2;
    shoulderCenterDrop = [0, 0, -S2sc];
    upperback = (mean(trc.data.SL(crop1:crop2,:)) + mean(trc.data.SR(crop1:crop2,:)))/2 + shoulderCenterDrop;
    lshoulder = mean(trc.data.SL(crop1:crop2,:)) + shoulderCenterDrop;
    rshoulder = mean(trc.data.SR(crop1:crop2,:)) + shoulderCenterDrop;
    lelbow = (mean(trc.data.ELLat(crop1:crop2,:)) + mean(trc.data.ELMed(crop1:crop2,:)))/2;
    relbow = (mean(trc.data.ERLat(crop1:crop2,:)) + mean(trc.data.ERMed(crop1:crop2,:)))/2;
    lwrist = (mean(trc.data.WLLat(crop1:crop2,:)) + mean(trc.data.WLMed(crop1:crop2,:)))/2;
    rwrist = (mean(trc.data.WRLat(crop1:crop2,:)) + mean(trc.data.WRMed(crop1:crop2,:)))/2;
    
    
%% Create model file and attach markers
%     addpath(cd(cd('..')));
    addpath(EKFCodePath);
    
    model = rlLieModel('JumpModel_Lie.xml');
    transformNames = {model.transforms.name};
    model.forwardPosition;
    
    initPos = zeros(size(model.position));
    % Joint position labels: 
    % (1-6) p0,p1,p2,r0,r1,r2, 
    % (7-9) backFB, backAxial, backLateral, 
    % (10-14) rshldrElev, rshldrAbd, rshldrExtRot, relbowFlex, relbowSup
    % (15-19) lshldrElev, lshldrAbd, lshldrExtRot, lelbowFlex, lelbowSup
    % (20-22) rhipFlex, rhipAbd, rhipExtRot, 
    % (23-26) rkneeExtend, rkneeExtRot, rankleDorsi, ranklePron
    % (27-29) lhipFlex, lhipAbd, lhipExtRot, 
    % (30-33) lkneeExtend, lkneeExtRot, lankleDorsi, lanklePron
    
    %Calculate init joint positions of shoulder abd. and forearm supination
    %(assume torso facing in global +X direction)
    initPos(11) = atan2(-(relbow(2)-rshoulder(2)),-(relbow(3)-rshoulder(3))); %shoulder abduction
    initPos(12) = -pi/4; %45 deg shoulder internal rotation, since ppl have tendency to do this during T-pose
    initPos(13) = deg2rad(5); %put small bend in elbow in correct direction
    initPos(14) = -pi/2; %assume forearm supinated 90 degrees (palm facing down in T-pose)
    initPos(16) = atan2(lelbow(2)-lshoulder(2),-(lelbow(3)-lshoulder(3))); %shoulder abduction
    initPos(17) = -pi/4; %45 deg shoulder internal rotation, since ppl have tendency to do this during T-pose
    initPos(18) = deg2rad(5); %put small bend in elbow in correct direction
    initPos(19) = -pi/2; %assume forearm supinated 90 degrees (palm facing down in T-pose)
    initPos(23) = deg2rad(-5); %put small bend in knee in correct direction
    initPos(30) = deg2rad(-5); %put small bend in knee in correct direction
    
     
    initPos_lie = model.position.copy;
    % Convert euler angle positions above to lie joint mapping
    R = angle2dcm(initPos(11),0,initPos(12),'XYZ');
    initPos_lie(3) = SO3(R'); %right shoulder
    initPos_lie(4) = SO2(initPos(13)); %right elbow
    initPos_lie(5) = SO2(initPos(14)); %right forearm
    R = angle2dcm(-initPos(16),0,initPos(17),'XYZ');
    initPos_lie(6) = SO3(R'); %left shoulder
    initPos_lie(7) = SO2(initPos(18)); %left elbow
    initPos_lie(8) = SO2(initPos(19)); %left forearm
    initPos_lie(10) = SO2(initPos(23)); %right knee
    initPos_lie(15) = SO2(initPos(30)); %left knee
    
    if(visu)
        vis = rlVisualizer('visModel',640,960);
        vis.addModel(model);
        vis.update();
    end
    
    %Save marker attachment frames and offsets
    markerFrames = {};
    markerOffsets = repmat(eye(4),1,1,30);
    
    %World to pelvis center
    trNum = find(ismember(transformNames,'world_to_base')==1);
    model.transforms(trNum).t = world_to_base;
    model_transforms.world_to_base = world_to_base;
    model.forwardPosition;
%     vis.update();
    
    % Attach ASIS and back markers
    trNum = find(ismember(transformNames,'midAsis2Back')==1); %last rotation joint to midAsis
    FT = model.transforms(trNum).frame_in.t;
    M6 = LieSensorCore('BT');
    M7 = LieSensorCore('BR');
    M8 = LieSensorCore('BL');
    markerOffsets(1:3,4,6) = FT(1:3,1:3)'*...
        (mean(trc.data.BT(crop1:crop2,:)) - FT(1:3,4)')';
    markerOffsets(1:3,4,7) = FT(1:3,1:3)'*...
        (mean(trc.data.BR(crop1:crop2,:)) - FT(1:3,4)')';
    markerOffsets(1:3,4,8) = FT(1:3,1:3)'*...
        (mean(trc.data.BL(crop1:crop2,:)) - FT(1:3,4)')';
    markerFrames{6} = model.transforms(trNum).frame_in.name; % note this is frame_IN, b/c want to attach to midAsis, not back0   
    markerFrames{7} = model.transforms(trNum).frame_in.name;
    markerFrames{8} = model.transforms(trNum).frame_in.name;
    
    
    %Back, shoulders, elbows, wrists
    trNum = find(ismember(transformNames,'back2Upperback')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t;
%     T(1:3,4) = FT(1:3,1:3)'*(mean(midShoulders(cropSta:cropEnd,:))-FT(1:3,4)')';
    T(1:3,4) = FT(1:3,1:3)'*([0,0,norm(upperback - FT(1:3,4)')])';
    model.transforms(trNum).t = T;
    model_transforms.back2Upperback = T;
    model.forwardPosition;
%     vis.update();
    M1 = LieSensorCore('Neck');
    FT2 = model.transforms(trNum).frame_out.t;
%     markerOffsets(1:3,4,1) = FT(1:3,1:3)'*rotx(-pi/2)*...
%         (mean(trc.data.Neck(crop1:crop2,:)) - FT2(1:3,4)')';
    markerOffsets(1:3,4,1) = FT2(1:3,1:3)'*...
        (mean(trc.data.Neck(crop1:crop2,:)) - FT2(1:3,4)')';
%     markerOffsets(1:3,4,1) = mean(trc.data.Neck(crop1:crop2,:))' - FT2(1:3,4);
    markerFrames{1} = model.transforms(trNum).frame_out.name;
    
    trNum = find(ismember(transformNames,'upperBack2rshoulder')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t;
    T(1:3,4) = FT(1:3,1:3)'*(rshoulder - FT(1:3,4)')';
    model.transforms(trNum).t = T;
    model_transforms.upperBack2rshoulder = T;
    model.forwardPosition;
%     vis.update();
    M2 = LieSensorCore('SR');
    FT2 = model.transforms(trNum).frame_out.t;
%     markerOffsets(1:3,4,2) = FT(1:3,1:3)'*rotx(-pi/2)*...
%         (mean(trc.data.SR(crop1:crop2,:)) - FT(1:3,4)')';
    markerOffsets(1:3,4,2) = FT2(1:3,1:3)'*...
        (mean(trc.data.SR(crop1:crop2,:)) - FT2(1:3,4)')';
%     markerOffsets(1:3,4,2) = mean(trc.data.SR(crop1:crop2,:))' - FT(1:3,4);
    markerFrames{2} = model.transforms(trNum).frame_out.name;
    
    trNum = find(ismember(transformNames,'upperBack2lshoulder')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t;
    T(1:3,4) = FT(1:3,1:3)'*(lshoulder - FT(1:3,4)')';
    model.transforms(trNum).t = T;
    model_transforms.upperBack2lshoulder = T;
    model.forwardPosition;
%     vis.update();
    M3 = LieSensorCore('SL');
    FT2 = model.transforms(trNum).frame_out.t;
%     markerOffsets(1:3,4,3) = FT(1:3,1:3)'*rotx(-pi/2)*...
%         (mean(trc.data.SL(crop1:crop2,:)) - FT(1:3,4)')';
    markerOffsets(1:3,4,3) = FT2(1:3,1:3)'*...
        (mean(trc.data.SL(crop1:crop2,:)) - FT2(1:3,4)')';
%     markerOffsets(1:3,4,3) = mean(trc.data.SL(crop1:crop2,:))' - FT(1:3,4);
    markerFrames{3} = model.transforms(trNum).frame_out.name;
    
    trNum = find(ismember(transformNames,'rShoulder2Elbow')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t;
%     T(1:3,4) = FT(1:3,1:3)'*(rshoulder - FT(1:3,4)')';
    T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(rshoulder - relbow)]';
    model.transforms(trNum).t = T;
    model_transforms.rShoulder2Elbow = T;
    model.forwardPosition;
%     vis.update();
    M9 = LieSensorCore('ERLat');
    M10 = LieSensorCore('ERMed');
    FT2 = model.transforms(trNum).frame_out.t;
%     markerOffsets(1:3,4,9) = FT2(1:3,1:3)'*...
%         (mean(trc.data.ERLat(crop1:crop2,:)) - FT2(1:3,4)')';
%     markerOffsets(1:3,4,10) = FT2(1:3,1:3)'*...
%         (mean(trc.data.ERMed(crop1:crop2,:)) - FT2(1:3,4)')';
%     markerDist = norm( mean(trc.data.ERLat(crop1:crop2,:)) - ...
%                        mean(trc.data.ERMed(crop1:crop2,:)))/2;
    markerDistLat = norm( mean(trc.data.ERLat(crop1:crop2,:)) - relbow);
    markerDistMed = norm( mean(trc.data.ERMed(crop1:crop2,:)) - relbow);
    if(markerDistLat<0.06 && markerDistMed>0.06) %for trc. data errors
        markerDistMed = markerDistLat;
    elseif(markerDistLat>0.06 && markerDistMed<0.06)
        markerDistLat = markerDistMed;
    elseif(markerDistLat>0.06 && markerDistMed>0.06)   
        markerDistLat = 0.035; %if both distances incorrect, assume this
        markerDistMed = 0.035;
    end
    markerOffsets(1:3,4,9) = FT2(1:3,1:3)'*[0, -markerDistLat, 0]';
    markerOffsets(1:3,4,10) = FT2(1:3,1:3)'*[0, markerDistMed, 0]';
    markerFrames{9} = model.transforms(trNum).frame_out.name;
    markerFrames{10} = model.transforms(trNum).frame_out.name;
    
    trNum = find(ismember(transformNames,'lShoulder2Elbow')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t;
%     T(1:3,4) = FT(1:3,1:3)'*(rshoulder - FT(1:3,4)')';
    T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(lshoulder - lelbow)]';
    model.transforms(trNum).t = T;
    model_transforms.lShoulder2Elbow = T;
    model.forwardPosition;
%     vis.update();
    M13 = LieSensorCore('ELLat');
    M14 = LieSensorCore('ELMed');
    FT2 = model.transforms(trNum).frame_out.t;
%     markerOffsets(1:3,4,13) = FT2(1:3,1:3)'*...
%         (mean(trc.data.ELLat(crop1:crop2,:)) - FT2(1:3,4)')';
%     markerOffsets(1:3,4,14) = FT2(1:3,1:3)'*...
%         (mean(trc.data.ELMed(crop1:crop2,:)) - FT2(1:3,4)')';
    % Above doesn't work because of inaccurate T-pose offset, assume elbow
    % and wrist markers equally on each side of joint center
    %     markerDist = norm( mean(trc.data.ERLat(crop1:crop2,:)) - ...
%                        mean(trc.data.ERMed(crop1:crop2,:)))/2;
    markerDistLat = norm( mean(trc.data.ELLat(crop1:crop2,:)) - lelbow);
    markerDistMed = norm( mean(trc.data.ELMed(crop1:crop2,:)) - lelbow);
    if(markerDistLat<0.06 && markerDistMed>0.06) %for trc. data errors
        markerDistMed = markerDistLat;
    elseif(markerDistLat>0.06 && markerDistMed<0.06)
        markerDistLat = markerDistMed;
    elseif(markerDistLat>0.06 && markerDistMed>0.06)   
        markerDistLat = 0.035; %if both distances incorrect, assume this
        markerDistMed = 0.035;
    end
    markerOffsets(1:3,4,13) = FT2(1:3,1:3)'*[0, markerDistLat, 0]';
    markerOffsets(1:3,4,14) = FT2(1:3,1:3)'*[0, -markerDistMed, 0]';
    markerFrames{13} = model.transforms(trNum).frame_out.name;
    markerFrames{14} = model.transforms(trNum).frame_out.name;
    
    trNum = find(ismember(transformNames,'rElbow2Wrist')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t;
%     T(1:3,4) = FT(1:3,1:3)'*(rshoulder - FT(1:3,4)')';
    T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(relbow - rwrist)]';
    model.transforms(trNum).t = T;
    model_transforms.rElbow2Wrist = T;
    model.forwardPosition;
%     vis.update();
    M11 = LieSensorCore('WRLat');
    M12 = LieSensorCore('WRMed');
    FT2 = model.transforms(trNum).frame_out.t;
%     markerOffsets(1:3,4,11) = FT2(1:3,1:3)'*...
%         (mean(trc.data.WRLat(crop1:crop2,:)) - FT2(1:3,4)')';
%     markerOffsets(1:3,4,12) = FT2(1:3,1:3)'*...
%         (mean(trc.data.WRMed(crop1:crop2,:)) - FT2(1:3,4)')';
% 	markerDist = norm( mean(trc.data.WRLat(crop1:crop2,:)) - ...
%                        mean(trc.data.WRMed(crop1:crop2,:)))/2;
    markerDistLat = norm( mean(trc.data.WRLat(crop1:crop2,:)) - rwrist);
    markerDistMed = norm( mean(trc.data.WRMed(crop1:crop2,:)) - rwrist);
    if(markerDistLat<0.06 && markerDistMed>0.06) %for trc. data errors
        markerDistMed = markerDistLat;
    elseif(markerDistLat>0.06 && markerDistMed<0.06)
        markerDistLat = markerDistMed;
    elseif(markerDistLat>0.06 && markerDistMed>0.06)   
        markerDistLat = 0.035; %if both distances incorrect, assume this
        markerDistMed = 0.035;
    end
    markerOffsets(1:3,4,11) = FT2(1:3,1:3)'*[0, -markerDistLat, 0]';
    markerOffsets(1:3,4,12) = FT2(1:3,1:3)'*[0, markerDistMed, 0]';
    markerFrames{11} = model.transforms(trNum).frame_out.name;
    markerFrames{12} = model.transforms(trNum).frame_out.name;
    
    trNum = find(ismember(transformNames,'lElbow2Wrist')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t;
%     T(1:3,4) = FT(1:3,1:3)'*(rshoulder - FT(1:3,4)')';
    T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(lelbow - lwrist)]';
    model.transforms(trNum).t = T;
    model_transforms.lElbow2Wrist = T;
    model.forwardPosition;
%     vis.update();
    M15 = LieSensorCore('WLLat');
    M16 = LieSensorCore('WLMed');
    FT2 = model.transforms(trNum).frame_out.t;
%     markerOffsets(1:3,4,15) = FT2(1:3,1:3)'*...
%         (mean(trc.data.WLLat(crop1:crop2,:)) - FT2(1:3,4)')';
%     markerOffsets(1:3,4,16) = FT2(1:3,1:3)'*...
%         (mean(trc.data.WLMed(crop1:crop2,:)) - FT2(1:3,4)')';
%     markerDist = norm( mean(trc.data.WLLat(crop1:crop2,:)) - ...
%                        mean(trc.data.WLMed(crop1:crop2,:)))/2;
    markerDistLat = norm( mean(trc.data.WLLat(crop1:crop2,:)) - lwrist);
    markerDistMed = norm( mean(trc.data.WLMed(crop1:crop2,:)) - lwrist);
    if(markerDistLat<0.06 && markerDistMed>0.06) %for trc. data errors
        markerDistMed = markerDistLat;
    elseif(markerDistLat>0.06 && markerDistMed<0.06)
        markerDistLat = markerDistMed;
    elseif(markerDistLat>0.06 && markerDistMed>0.06)   
        markerDistLat = 0.035; %if both distances incorrect, assume this
        markerDistMed = 0.035;
    end
    markerOffsets(1:3,4,15) = FT2(1:3,1:3)'*[0, markerDistLat, 0]';
    markerOffsets(1:3,4,16) = FT2(1:3,1:3)'*[0, -markerDistMed, 0]';
    markerFrames{15} = model.transforms(trNum).frame_out.name;
    markerFrames{16} = model.transforms(trNum).frame_out.name;
    
    
    %Hips, knees, ankles, feet
    trNum = find(ismember(transformNames,'rAsis2Hip')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t; %mid_asis
    T(1:3,4) = mid2RHC;
    model.transforms(trNum).t = T;
    model_transforms.rAsis2Hip = T;
    model.forwardPosition;
%     vis.update();
    M4 = LieSensorCore('HR');
    markerOffsets(1:3,4,4) = FT(1:3,1:3)'*...
        (mean(trc.data.HR(crop1:crop2,:)) - FT(1:3,4)')';
%     markerOffsets(1:3,4,4) = mean(trc.data.HR(crop1:crop2,:))' - FT(1:3,4);
    markerFrames{4} = model.transforms(trNum).frame_in.name;
    
    trNum = find(ismember(transformNames,'lAsis2Hip')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t; %mid_asis
    T(1:3,4) = mid2LHC;
    model.transforms(trNum).t = T;
    model_transforms.lAsis2Hip = T;
    model.forwardPosition;
%     vis.update();
    M5 = LieSensorCore('HL');
    markerOffsets(1:3,4,5) = FT(1:3,1:3)'*...
        (mean(trc.data.HL(crop1:crop2,:)) - FT(1:3,4)')';
%     markerOffsets(1:3,4,5) = mean(trc.data.HL(crop1:crop2,:))' - FT(1:3,4);
    markerFrames{5} = model.transforms(trNum).frame_in.name;
    
    trNum = find(ismember(transformNames,'rHip2Knee')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t;
    T(1:3,4) = FT(1:3,1:3)'*(rknee - FT(1:3,4)')';
%     T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm((mid+mid2RHC) - rknee)]';
    model.transforms(trNum).t = T;
    model_transforms.rHip2Knee = T;
    model.forwardPosition;
%     vis.update();
    M17 = LieSensorCore('KRLat');
    M18 = LieSensorCore('KRMed');
    FT2 = model.transforms(trNum).frame_out.t;
    markerOffsets(1:3,4,17) = FT2(1:3,1:3)'*...
        (mean(trc.data.KRLat(crop1:crop2,:)) - rknee)';
    markerOffsets(1:3,4,18) = FT2(1:3,1:3)'*...
        (mean(trc.data.KRMed(crop1:crop2,:)) - rknee)';
%     markerFrames{17} = model.transforms(trNum).frame_in.name;
%     markerFrames{18} = model.transforms(trNum).frame_in.name;
    markerFrames{17} = model.transforms(trNum).frame_out.name;
    markerFrames{18} = model.transforms(trNum).frame_out.name;
    
    trNum = find(ismember(transformNames,'lHip2Knee')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t;
    T(1:3,4) = FT(1:3,1:3)'*(lknee - FT(1:3,4)')';
%     T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm((mid+mid2LHC) - lknee)]';
    model.transforms(trNum).t = T;
    model_transforms.lHip2Knee = T;
    model.forwardPosition;
%     vis.update();
    M24 = LieSensorCore('KLLat');
    M25 = LieSensorCore('KLMed');
    FT2 = model.transforms(trNum).frame_out.t;
    markerOffsets(1:3,4,24) = FT2(1:3,1:3)'*...
        (mean(trc.data.KLLat(crop1:crop2,:)) - lknee)';
    markerOffsets(1:3,4,25) = FT2(1:3,1:3)'*...
        (mean(trc.data.KLMed(crop1:crop2,:)) - lknee)';
%     markerFrames{24} = model.transforms(trNum).frame_in.name;
%     markerFrames{25} = model.transforms(trNum).frame_in.name;
    markerFrames{24} = model.transforms(trNum).frame_out.name;
    markerFrames{25} = model.transforms(trNum).frame_out.name;
    
    trNum = find(ismember(transformNames,'rKnee2Ankle')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t;
%     T(1:3,4) = FT(1:3,1:3)'*(mean(rankle(crop1:crop2,:))-FT(1:3,4)')';
    T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(rknee - rankle)]';
    model.transforms(trNum).t = T;
    model_transforms.rKnee2Ankle = T;
    model.forwardPosition;
%     vis.update();
    M19 = LieSensorCore('ARLat');
    M20 = LieSensorCore('ARMed');
    FT2 = model.transforms(trNum).frame_out.t;
    markerOffsets(1:3,4,19) = FT2(1:3,1:3)'*...
        (mean(trc.data.ARLat(crop1:crop2,:)) - rankle)';
    markerOffsets(1:3,4,20) = FT2(1:3,1:3)'*...
        (mean(trc.data.ARMed(crop1:crop2,:)) - rankle)';
%     markerFrames{19} = model.transforms(trNum).frame_in.name;
%     markerFrames{20} = model.transforms(trNum).frame_in.name;
    markerFrames{19} = model.transforms(trNum).frame_out.name;
    markerFrames{20} = model.transforms(trNum).frame_out.name;
    
    trNum = find(ismember(transformNames,'lKnee2Ankle')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t;
%     T(1:3,4) = FT(1:3,1:3)'*(mean(rankle(crop1:crop2,:))-FT(1:3,4)')';
    T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(lknee - lankle)]';
    model.transforms(trNum).t = T;
    model_transforms.lKnee2Ankle = T;
    model.forwardPosition;
%     vis.update();
    M26 = LieSensorCore('ALLat');
    M27 = LieSensorCore('ALMed');
    FT2 = model.transforms(trNum).frame_out.t;
    markerOffsets(1:3,4,26) = FT2(1:3,1:3)'*...
        (mean(trc.data.ALLat(crop1:crop2,:)) - lankle)';
    markerOffsets(1:3,4,27) = FT2(1:3,1:3)'*...
        (mean(trc.data.ALMed(crop1:crop2,:)) - lankle)';
%     markerFrames{26} = model.transforms(trNum).frame_in.name;
%     markerFrames{27} = model.transforms(trNum).frame_in.name;
    markerFrames{26} = model.transforms(trNum).frame_out.name;
    markerFrames{27} = model.transforms(trNum).frame_out.name;
    
    trNum = find(ismember(transformNames,'rAnkle2Toe')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t;
    T(1:3,4) = FT(1:3,1:3)'*(rtoe - FT(1:3,4)')';
%     T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(rankle - rtoe)]';
    model.transforms(trNum).t = T;
    model_transforms.rAnkle2Toe = T;
    model.forwardPosition;
%     vis.update();
    M21 = LieSensorCore('FRLat');
    M22 = LieSensorCore('FRMed');
    M23 = LieSensorCore('FRHeel');
    FT2 = model.transforms(trNum).frame_out.t;
    markerOffsets(1:3,4,21) = FT2(1:3,1:3)'*...
        (mean(trc.data.FRLat(crop1:crop2,:)) - rtoe)';
    markerOffsets(1:3,4,22) = FT2(1:3,1:3)'*...
        (mean(trc.data.FRMed(crop1:crop2,:)) - rtoe)';
    markerOffsets(1:3,4,23) = FT2(1:3,1:3)'*...
        (mean(trc.data.FRHeel(crop1:crop2,:)) - rtoe)';
%     markerFrames{21} = model.transforms(trNum).frame_in.name;
%     markerFrames{22} = model.transforms(trNum).frame_in.name;
%     markerFrames{23} = model.transforms(trNum).frame_in.name;
    markerFrames{21} = model.transforms(trNum).frame_out.name;
    markerFrames{22} = model.transforms(trNum).frame_out.name;
    markerFrames{23} = model.transforms(trNum).frame_out.name;
    
    trNum = find(ismember(transformNames,'lAnkle2Toe')==1);
    T = model.transforms(trNum).t;
    FT = model.transforms(trNum).frame_in.t;
    T(1:3,4) = FT(1:3,1:3)'*(ltoe - FT(1:3,4)')';
%     T(1:3,4) = FT(1:3,1:3)'*[0,0,-norm(lankle - ltoe)]';
    model.transforms(trNum).t = T;
    model_transforms.lAnkle2Toe = T;
    model.forwardPosition;
%     vis.update();
    M28 = LieSensorCore('FLLat');
    M29 = LieSensorCore('FLMed');
    M30 = LieSensorCore('FLHeel');
    FT2 = model.transforms(trNum).frame_out.t;
    markerOffsets(1:3,4,28) = FT2(1:3,1:3)'*...
        (mean(trc.data.FLLat(crop1:crop2,:)) - ltoe)';
    markerOffsets(1:3,4,29) = FT2(1:3,1:3)'*...
        (mean(trc.data.FLMed(crop1:crop2,:)) - ltoe)';
    markerOffsets(1:3,4,30) = FT2(1:3,1:3)'*...
        (mean(trc.data.FLHeel(crop1:crop2,:)) - ltoe)';
%     markerFrames{28} = model.transforms(trNum).frame_in.name;
%     markerFrames{29} = model.transforms(trNum).frame_in.name;
%     markerFrames{30} = model.transforms(trNum).frame_in.name;
    markerFrames{28} = model.transforms(trNum).frame_out.name;
    markerFrames{29} = model.transforms(trNum).frame_out.name;
    markerFrames{30} = model.transforms(trNum).frame_out.name;
    
    
    if(visu)
        vis.update();
        clear vis;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Attach markers to model
    % Order: HR, HL, R_shldr, L_shldr, knee_lat, knee_med,
    % ankle_lat, ankle_med
    
    
    if(attachMarkers)
        Markers = [M1 M2 M3 M4 M5 M6 M7 M8 M9 M10 M11 M12 M13 M14 M15...
                   M16 M17 M18 M19 M20 M21 M22 M23 M24 M25 M26 M27 M28 M29 M30];
        
        for i = 1:numel(Markers)
            Markers(i).addDecorator('position');
            model.addSensor(Markers(i), markerFrames{i}, markerOffsets(:,:,i));
            
%             if(visu)
%                 model.forwardPosition();
%                 vis = rlVisualizer('visModel',640,960);
%                 vis.addModel(model);
%                 vis.update();
%                 clear vis;
%             end
        end
        model.forwardPosition();
        
        if(visu)
            model.forwardPosition();
            vis = rlVisualizer('visModel',640,960);
            vis.addModel(model);
            vis.update();
            clear vis;
        end
    end
    
end

