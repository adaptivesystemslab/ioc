% clear in correct order
% close all;
clear vis;
clear model;
clearvars;

visualize = 1;

filepathDll = '.';
filepathModel = fullfile('.', 'Models', '7dof_squat.xml');
pathToRawData = fullfile('D:', 'aslab', 'data', 'Squats_TUAT_2015-12', 'Subject01', 'Session1', 'SQUA_STD_NON2');
addpath(fullfile('..', 'Common'));
cd(filepathDll);

% load the kinematic data	
filesToLoad{1} = fullfile(pathToRawData, searchForFileByExt(pathToRawData, 'Kinematics_meas*.mat'));
filesToLoad{2} = fullfile(pathToRawData, searchForFileByExt(pathToRawData, 'Dynamics_meas_*.mat'));

load1 = load(filesToLoad{1}); % q, dq, ddq, Markers
load3 = load(filesToLoad{2}); % body_height, body_weight, FP, COP

dt = 0.01;
fullDataAngles = load1.q;

% q = [0 0 pi/2 -pi/2 pi/2 -pi/2 0]';
% q = [0 0 0 0 0 0 0]';
qIn = fullDataAngles;
% q = fullDataAngles(:, 1);
dqIn = calcDeriv(qIn, dt);
ddqIn = calcDeriv(dqIn, dt);
% dddq = calcDeriv(ddq, dt);

% % q = [0 pi/2 pi/8 pi/8 pi/8 pi/8 pi/8]';

% create the model and set it to default position
model = rlCModel(filepathModel);
model.forwardPosition();

%Add accelerometer to the end efffector
sens = SensorCore('sens');
sens.addDecorator('accelerometer');
model.addSensor(sens,'frame8',eye(4));

% update linkage data
% param.body_weight = 60.9;
% param.body_height = 1.68;
param.body_height=1.75;%Sb 6
param.body_weight=75.3;
    
param = updateLinkInfo(param, load1.Markers);
model = updateModelInfo(model, param);

model.position(:) = zeros(size(qIn(:, 1)));
model.forwardPosition();

% model.bodies(:).t

% create visualizer window and add the model to it
if visualize
%     vis = rlVisualizer('vis',640,480);
    vis = rlVisualizer('vis',640*2,480*2);
    vis.addModel(model);
    vis.update();
end

q = [];
dq = [];
ddq = [];
dddq = [];
x = [];
dx = [];
ddx = [];
dddx = [];
tau = [];
dtau = [];
ddtau = [];
cop = [];
dcop = [];
ddcop = [];
com = [];
dcom = [];
ddcom = [];
ep = [];
ek = [];
geo = [];
en = [];

efToUse = 5;

for i=1:size(qIn, 2)
    if mod(i, 100) == 0
        fprintf('Frame: %u/%u\n', i, size(qIn, 2));
    end
    
    %SET JOINT ANGLES HERE
    
    %THIS IS TO SET THEM IN ORDER
    model = inputData(model, qIn(:, i)', dqIn(:, i)', ddqIn(:, i)');
    
    %MUST CALL THIS TO UPDATE MODEL
    model.forwardPosition();
    model.forwardVelocity();
    model.forwardAcceleration();
    model.inverseDynamics();
    
    f = model.getFrameByName('frame8');
    ddx_curr = f.a(4:6) + cross(f.v(1:3),f.v(4:6)) - f.t(1:3,1:3)'*model.g;
    ddx_curr = rotz(pi/2)*ddx_curr;
    tau_curr = model.torque;
    
    model.calculateJacobian();              model = inputData(model, qIn(:, i)', dqIn(:, i)', ddqIn(:, i)');
    model.calculateJacobianDerivative();    model = inputData(model, qIn(:, i)', dqIn(:, i)', ddqIn(:, i)');
    
    x_curr = model.bodies(efToUse).t(1:3, 4);
    dx_curr = model.J*model.velocity;
    dx_curr = dx_curr([1 2 3]);
    
        
 
%     ddx_curr = model.J*model.acceleration + model.JdQd;
%     ddx1 = model.J*model.acceleration;
%     ddx2 = model.JdQd;
%     ddx_curr = ddx_curr([1 2 3]);

%     q_curr = model.position(:);
%     dq_curr = model.velocity(:);
%     ddq_curr = model.acceleration(:);
   
%     model.calculateMassMatrix();            model = inputData(model, qIn(:, i)', dqIn(:, i)', ddqIn(:, i)');
%     model.calculateCentrifugalCoriolis();   model = inputData(model, qIn(:, i)', dqIn(:, i)', ddqIn(:, i)');
%     model.calculateGravity();               model = inputData(model, qIn(:, i)', dqIn(:, i)', ddqIn(:, i)');


    
%     M(:, :, i) = model.M;
%     q =     [q      q_curr];
%     dq =    [dq     dq_curr];
%     ddq =   [ddq    ddq_curr];
    x =     [x      x_curr];
    dx =    [dx     dx_curr];
    ddx =   [ddx    ddx_curr];
    tau =   [tau    tau_curr];
    
    if visualize
        %UPDATE VIS AND CAPTURE FRAME
        vis.update();
        
        pause(0.01);
    end
end


% figure;hold on;
% plot(q(1,:),'r')
% plot(q(2,:),'g')
% plot(q(3,:),'b')
% plot(q(4,:),'y')
% plot(q(5,:),'m')
% plot(q(6,:),'k')
% plot(q(7,:),'c')
% ylabel('Joint angle [rad]')
% 
% figure;subplot(211);hold on;
% plot(dq(1,:),'r')
% plot(dq(2,:),'g')
% plot(dq(3,:),'b')
% plot(dq(4,:),'y')
% plot(dq(5,:),'m')
% plot(dq(6,:),'k')
% plot(dq(7,:),'c')
% ylabel('Joint velocity [rad.s-1]')
% subplot(212);hold on;
% plot(ddq(1,:),'r')
% plot(ddq(2,:),'g')
% plot(ddq(3,:),'b')
% plot(ddq(4,:),'y')
% plot(ddq(5,:),'m')
% plot(ddq(6,:),'k')
% plot(ddq(7,:),'c')
% ylabel('Joint acc [rad.s-2]')

figure
subplot(211);hold on;
plot(x(1,:),'r')
plot(x(2,:),'g')
plot(x(3,:),'b')
ylabel('Linear pos')
subplot(212);hold on;
plot(dx(1,:),'r')
plot(dx(2,:),'g')
plot(dx(3,:),'b')
ylabel('Linear vel')

figure
h1 = subplot(211);hold on;
plot(ddx(1,:),'r')
plot(ddx(2,:),'g')
plot(ddx(3,:),'b')
ylabel('Linear acc calc')
h2 = subplot(212);hold on;
ddx2 = calcDeriv(dx, 0.01);
plot(ddx2(1,:),'r')
plot(ddx2(2,:),'g')
plot(ddx2(3,:),'b')
ylabel('Linear acc diff')
linkaxes([h1 h2], 'xy');

% figure;hold on;
% plot(tau(1,:),'r')
% plot(tau(2,:),'g')
% plot(tau(3,:),'b')
% plot(tau(4,:),'y')
% plot(tau(5,:),'m')
% plot(tau(6,:),'k')
% ylabel('Joint Torques [N.m]')
% title('RL formulation');

% figure;
% plot(endeffpos(:, (7:9)+6)); title('knee end eff');
