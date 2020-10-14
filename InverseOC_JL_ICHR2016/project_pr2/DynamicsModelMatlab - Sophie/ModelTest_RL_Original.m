% clear vis
% % clear mdl
clc

addpath(genpath('DynamicsModelMatlab'));

%% General usage of RL wrapper
% load model
mdl = rlCModel('robot.rlmdl.xml');
% apply forward kinematics to determine joint positions 
mdl.forwardPosition();

% for all joints, print their position
for i = 1:length(mdl.position)
    fprintf('Position of joint %s is [%.2f, %.2f, %.2f]\n', mdl.joints(i).name, mdl.joints(i).frame_in.t(1:3, 4));
end

fprintf('\n')

% Apply a different joint configuration
mdl.position = [pi/4, 0, pi/4, 0, 0, 0, 0];
mdl.forwardPosition();
for i = 1:length(mdl.position)
    fprintf('Position of joint %s is [%.2f, %.2f, %.2f]\n', mdl.joints(i).name, mdl.joints(i).frame_in.t(1:3, 4));
end

%% Example of how to use dynamic equations
% Generate a sequence of joint angles, velocities and accelerations
% Compute torque using inverse dynamics
%% Change q
dt = pi/100;
t = 0:dt:(6*pi);
w = 1;
qtraj = [sin(w*t); sin((1/2)*w*t); sin((3/2)*w*t); sin((1/2)*w*t); sin((3/2)*w*t); sin((1/2)*w*t); sin((3/2)*w*t)]';
dqtraj = calcDerivVert(qtraj, dt);
ddqtraj = calcDerivVert(dqtraj, dt);

for i = 1:length(t)
    q = qtraj(i, :);
    dq = dqtraj(i, :);
    ddq = ddqtraj(i, :);
    
    mdl.position = q;
    mdl.velocity = dq;
    mdl.acceleration = ddq;
    mdl.forwardPosition();
    mdl.forwardVelocity();
    mdl.forwardAcceleration();
    mdl.inverseDynamics();
    tau_rl(i, :) = mdl.torque';
end

plot(tau_rl)


%% Here is an example of how to use the visualizer

% vis = rlVisualizer('vis',640,480);
% vis.addModel(mdl);
% 
% colours.axis                 = [0 0 0 0.2];    % axis labels
% vis.addMarker('x-axis', [1 0 0], colours.axis);
% vis.addMarker('y-axis', [0 1 0], colours.axis);
% vis.addMarker('z-axis', [0 0 1], colours.axis);
% 
% vis.update();

%%
% mdl.position = zeros(size(mdl.position));
% for i = 1:length(mdl.position)
%     fprintf('Testing joint %u of %u: %s\n', i, length(mdl.position), mdl.joints(i).name);
%     
%     mdl.position = zeros(size(mdl.position));
%     q = [0:0.001:pi/2 pi/2:-0.001:0];
%     for j = 1:length(q)
%         mdl.position(i) = q(j);
%         mdl.forwardPosition();
%         vis.update();
%     end
% end