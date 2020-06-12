% clear vis
% % clear mdl
clc

% test ars_fullbody
addpath(genpath('D:\aslab_gitlab\kalmanfilter\ik_framework\ekf'));
t = eye(4);

mdl = rlCModel('.\model\ioc_v3_singleArm.xml');
mdl.forwardPosition();

vis = rlVisualizer('vis',640,480);
vis.addModel(mdl);

colours.axis                 = [0 0 0 0.2];    % axis labels
vis.addMarker('x-axis', [1 0 0], colours.axis);
vis.addMarker('y-axis', [0 1 0], colours.axis);
vis.addMarker('z-axis', [0 0 1], colours.axis);

vis.update();

%%
mdl.position = zeros(size(mdl.position));
for i = 1:length(mdl.position)
% for i = 6
    fprintf('Testing joint %u of %u: %s\n', i, length(mdl.position), mdl.joints(i).name);
    
    mdl.position = zeros(size(mdl.position));
    q = [0:0.001:pi/2 pi/2:-0.001:0];
    for j = 1:length(q)
        mdl.position(i) = q(j);
        mdl.forwardPosition();
        vis.update();
    end
end