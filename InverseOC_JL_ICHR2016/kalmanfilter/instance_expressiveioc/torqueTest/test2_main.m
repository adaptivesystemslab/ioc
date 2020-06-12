%% Set up Corke
cd('D:\aslab_gitlab\kalmanfilter\ik_framework\toolboxes\Robotics_Corke');
startup_rvc;

cd('D:\aslab_gitlab\kalmanfilter\ik_framework\instance_expressiveioc');
addpath(genpath('.\support'));
addpath(genpath('..\common'));
addpath(genpath('..\toolboxes'));
addpath(genpath('..\ekf'));

%% Set global parameters
g = -9.81;
m1 = 1.2540;
l1 = [0.3240  0       0];
% r1 = [0.1471 -0.0237 -0.0091];
r1 = [0.1471 0 0];
i1 = [    0.0038   -0.0001   -0.0026
   -0.0001    0.0143   -0.0003
   -0.0026   -0.0003    0.0143];

m2 = 0.7410;
l2 = [0.2192 0      0];
% r2 = [0.0901 0.0046 0.0042];
r2 = [0.0901 0 0];
i2 = [    0.0007   -0.0004    0.0006
   -0.0004    0.0024    0.0001
    0.0006    0.0001    0.0022];

%% Load Corke Model
l1c = rotz(pi)*l1';
r1c = rotz(-pi/2)*r1';
i1c = rotz(-pi/2)*i1*rotz(-pi/2)';
l2c = rotz(pi)*l2';
r2c = rotz(-pi/2)*r2';
i2c = rotz(-pi/2)*i2*rotz(-pi/2)';

L(1) = RevoluteMDH('d', 0,     'a', 0, 'alpha',     0, 'offset',    0, 'm', 0,  'r', [0 0 0],  'I', zeros(3,3),  'qlim', [-pi*2/3 pi/2]);
L(2) = RevoluteMDH('d', 0,     'a', 0, 'alpha', -pi/2, 'offset', pi/2, 'm', m1, 'r', r1c ,     'I', i1c,         'qlim', [-pi/2 pi/2]); 
L(3) = RevoluteMDH('d', l1(1), 'a', 0, 'alpha', pi/2,  'offset',    0, 'm', 0,  'r', [0 0 0],  'I', zeros(3,3),  'qlim', [-2*pi/3 pi]); 

L(4) = RevoluteMDH('d', 0,     'a', 0, 'alpha', -pi/2, 'offset', 0,    'm', m2, 'r', r2c,      'I', zeros(3, 3), 'qlim', [0 5*pi/6]); 
L(5) = RevoluteMDH('d', l2(1), 'a', 0, 'alpha', pi/2,  'offset', 0,    'm', 0,  'r', [0 0 0],  'I', i2,          'qlim', [0 0]);


mdl_corke = SerialLink(L, 'name', 'mdl_corke');

%% Load RL Model
filepathModelXml = fullfile('.', 'model', 'ioc_v3_dhAlign.xml');
mdl_rl = rlCModel(filepathModelXml);

allTransformNames = {mdl_rl.transforms.name};
allBodyNames = {mdl_rl.bodies.name};
 
currTransformStr = 'length_rshoulder_relbow';
currBodyStr = 'body_rshoulder_relbow';
indFrame = find(ismember(allTransformNames, currTransformStr) == 1);
T = eye(4);
T(1:3, 4) = l1;
mdl_rl.transforms(indFrame).t = T;

indFrame = find(ismember(allBodyNames, currBodyStr) == 1);
mdl_rl.bodies(indFrame).m   = m1;
mdl_rl.bodies(indFrame).com = r1;
mdl_rl.bodies(indFrame).I   = i1;

currTransformStr = 'length_relbow_rwrist';
currBodyStr = 'body_relbow_rwrist';
indFrame = find(ismember(allTransformNames, currTransformStr) == 1);
T = eye(4);
T(1:3, 1:3) = roty(pi/2);
T(1:3, 4) = l2;
mdl_rl.transforms(indFrame).t = T;

indFrame = find(ismember(allBodyNames, currBodyStr) == 1);
mdl_rl.bodies(indFrame).m   = m2;
mdl_rl.bodies(indFrame).com = r2;
mdl_rl.bodies(indFrame).I   = i2;
                
mdl_rl.forwardPosition();

%% Plot Models
vis = rlVisualizer('vis',640,480);
mdl_rl.forwardPosition();
vis.addModel(mdl_rl);
applyMarkersToVisualization(vis, mdl_rl, [], [], [], []);
vis.update();

%% Change q
dt = pi/100;
t = 0:dt:(6*pi);
w = 1;
qtraj = [sin(w*t); sin((1/2)*w*t); sin((3/2)*w*t); sin((1/2)*w*t); sin((3/2)*w*t)]';
dqtraj = calcDerivVert(qtraj, dt);
ddqtraj = calcDerivVert(dqtraj, dt);

% q = [0 0 0 0 0];
% dq = [0 0 0 0 0];
% ddq = [0 0 0 0 0];

for i = 1:length(t)
    q = qtraj(i, :);
    dq = dqtraj(i, :);
    ddq = ddqtraj(i, :);
    
    mdl_rl.position = q;
    mdl_rl.velocity = dq;
    mdl_rl.acceleration = ddq;
    mdl_rl.forwardPosition();
    mdl_rl.forwardVelocity();
    mdl_rl.forwardAcceleration();
    mdl_rl.inverseDynamics();
    
%     mdl_corke.plot(q);
%     vis.update();
    
    tau_corke(i, :) = mdl_corke.rne(q, dq, ddq);
    tau_rl(i, :) = mdl_rl.torque';
    tau_hand(i, :) = [0 m1*g*r1(1)+m2*g*(l1(1)+r2(1)) 0 m2*g*r2(1) 0];
    % tau_corke_rl = [tau_corke; tau_rl]
end
