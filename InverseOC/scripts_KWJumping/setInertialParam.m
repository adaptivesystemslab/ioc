function model = setInertialParam(model, trNum, gender, weight)


linkName = model.transforms(trNum).name;
linkLength = model.transforms(trNum).t(1:3,4);

% Set link type, apply rotations as required
% rotMtxDumas = rotz(pi/2)*rotx(pi/2);
rotMtxDumas = rotx(pi/2);
rotMirror = [1 0 0; 0 -1 0; 0 0 1];

switch linkName
    case {'upperBack2lshoulder','upperBack2rshoulder','lAsis2Hip'} % either doesn't exist or already accounted for in other limbs
        dumasFrameStr = '';

    case {'rAsis2Hip'}
        bodyName = 'mid_asis';
        dumasFrameStr = 'pelvis';
%         rotMtxDumas = [];

    case {'lHip2Knee', 'rHip2Knee'}
        bodyName = [linkName(1) 'hip5'];
        dumasFrameStr = 'thigh';

    case {'lKnee2Ankle', 'rKnee2Ankle'}
        bodyName = [linkName(1) 'knee5'];
        dumasFrameStr = 'leg';
%         rotMtxDumas = rotx(pi/2);

    case {'lAnkle2Toe', 'rAnkle2Toe'}
        bodyName = [linkName(1) 'ankle3'];
        dumasFrameStr = 'foot';
%         rotMtxDumas = [];

    case {'back2Upperback'}
        bodyName = ['back5'];
        dumasFrameStr = 'torso';
%         rotMtxDumas = roty(pi)*rotz(pi/2)*rotx(pi/2);
        rotMtxDumas = rotx(-pi/2);

%     case {'length_t1c7_c1head'}
%         dumasFrameStr = 'head&neck';

    case {'lShoulder2Elbow', 'rShoulder2Elbow'}
        bodyName = [linkName(1) 'shoulder5'];
        dumasFrameStr = 'arm';
%         rotMtxDumas = [];

    case {'lElbow2Wrist', 'rElbow2Wrist'}
        bodyName = [linkName(1) 'elbow3'];
        dumasFrameStr = 'forearm';
%         rotMtxDumas = [];

%     case {'length_rwrist_rhand', 'length_lwrist_lhand'}
%         dumasFrameStr = 'hand';
end


if(strcmpi(linkName(1),'l'))
    rotMtxDumas = rotMirror*rotMtxDumas;
end


mass =            lookupTableDumas('mass',dumasFrameStr, gender, [], [])*double(weight);
% Add extra mass for head and hands
% if(contains(linkName,'Wrist'))
%     addMass = lookupTableDumas('mass','hand', gender, [], [])*weight;
%     mass = mass + addMass;
% elseif(contains(linkName,'back2Upperback'))
%     addMass = lookupTableDumas('mass','head', gender, [], [])*weight;
%     mass = mass + addMass;
% end

comScale =        lookupTableDumas('com',dumasFrameStr, gender, norm(linkLength), []);
inertialScale =   lookupTableDumas('inertial',dumasFrameStr, gender, norm(linkLength), mass);

% parallel axis, since Dumas inertia assumes it's on the edge of
% the limb
% inertial_Hyugens =  mass * ...
%     [comScale(2)^2+comScale(3)^2,  -comScale(1)*comScale(2),     -comScale(1)*comScale(3);
%     -comScale(1)*comScale(2),       comScale(1)^2+comScale(3)^2, -comScale(2)*comScale(3);
%     -comScale(1)*comScale(3),      -comScale(2)*comScale(3),      comScale(1)^2+comScale(2)^2];

com = rotMtxDumas*comScale;
inertial = rotMtxDumas*inertialScale; %+ rotMtxDumas*inertial_Hyugens;

%Set body parameters

allBodyNames = {model.bodies.name};
bodyNum = find(ismember(allBodyNames, bodyName) == 1);

model.bodies(bodyNum).m = mass;
model.bodies(bodyNum).com = com;
model.bodies(bodyNum).I = inertial;