% RL setup testing

if 0
    % load the existing IIT model via getModel
    trialInfo.model = 'IIT_3DOF';
    trialInfo.baseModel = 'IIT';
    trialInfo.path = '..\Data\IK\2019_04_10_IIT_Squat02\matEkfIk\SQUA_STD_FAT_Subject01_mocap_mocap_X00_floating_ekfId1_ekfIk.mat';
    
    modelObj = getModel(trialInfo);
    traj = modelObj.loadData(trialInfo);
    model = modelObj.model;
    q = traj.q;
else
    % load the template class directly to verify implemented template is
    % working
    modelObj = ModelRL_Template();
    modelObj.loadModel('..\Libraries\rl\ik_framework\instance_iit\model\iit_v10.xml');
    model = modelObj.model;
    
    qIndiv = sin(0:0.01:8*pi);
    q = zeros(length(qIndiv), length(model.joints));
    q(:, 7) = qIndiv;
end

% visualize the loaded model
vis = rlVisualizer('vis',640,480);
model.forwardPosition();
vis.addModel(model);

vis.addMarker('x-axis', [1 0 0], [0.2 0.2 0.2 1]);
vis.addMarker('y-axis', [0 1 0], [0.2 0.2 0.2 1]);
vis.addMarker('z-axis', [0 0 1], [0.2 0.2 0.2 1]);

vis.update();

for i = 1:size(q, 1)
    model.position = q(i, :);
    model.forwardPosition();
    
    pause(0.001);
    vis.update();
end