%Test derivatives of transforms
addpath('..\');
%Create Model and visualize
mdl = rlCModel('simple_3dof_model.xml');
sens = SensorCore('sens1');
sens.addDecorator('gyroscope');
sens.addDecorator('accelerometer');

mdl.addSensor(sens,'frame_ee',eye(4));
mdl.forwardPosition();

vis = rlVisualizer('vis',640,480);
vis.addModel(mdl);
vis.update();

%% Numerically compute derivative of transforms
dTb_ee_num = zeros(4,4,numel(mdl.joints));
dTee_b_num = zeros(4,4,numel(mdl.joints));

ep = 1e-6;
for i=1:numel(mdl.joints)
    Tb_ee = sens.transform;
    Tee_b = SE3.fastInverse(Tb_ee);
    pos_init = mdl.position(i);
    mdl.position(i) = mdl.position(i) + ep;
    mdl.forwardPosition();
    Tb_ee_eps = sens.transform;
    Tee_b_eps = SE3.fastInverse(Tb_ee_eps);
    mdl.position(i)= pos_init;
    mdl.forwardPosition();
    
    dTb_ee_num(:,:,i) = (Tb_ee_eps-Tb_ee)/ep;
    dTee_b_num(:,:,i) = (Tee_b_eps-Tee_b)/ep;
    
end

%% Lets look at the derivative of one sensor in the other sensor frame

clear vis;
clear mdl;

mdl = rlCModel('..\Models\Lower_Body.xml');
mdl.forwardPosition();

%Add sensors on each ee
sens1 = SensorCore('sens1');
sens1.addDecorator('gyroscope');
sens1.addDecorator('accelerometer');
mdl.addSensor(sens1,'lankle3',eye(4));

sens2 = SensorCore('sens2');
sens2.addDecorator('gyroscope');
sens2.addDecorator('accelerometer');
mdl.addSensor(sens2,'rankle3',eye(4));
mdl.forwardPosition();

mdl.position = rand(numel(mdl.joints),1);
mdl.forwardPosition();

vis = rlVisualizer('vis',640,480);
vis.addModel(mdl);
vis.update();

dTee1_ee2_num = zeros(4,4,numel(mdl.joints));
ep = 1e-6;
for i=1:numel(mdl.joints)
    Tb_ee2 = sens2.transform;
    Tee1_b = SE3.fastInverse(sens1.transform);
    Tee1_ee2 = Tee1_b*Tb_ee2;
    pos_init = mdl.position(i);
    mdl.position(i) = mdl.position(i) + ep;
    mdl.forwardPosition();
    Tb_ee2_eps = sens2.transform;
    Tee1_b_eps = SE3.fastInverse(sens1.transform);
    Tee1_ee2_eps = Tee1_b_eps*Tb_ee2_eps;
    
    mdl.position(i)= pos_init;
    mdl.forwardPosition();
    dTee1_ee2_num(:,:,i) = (Tee1_ee2_eps-Tee1_ee2)/ep;
end

%% Now Analitically 

dTee1_ee2 = zeros(4,4,numel(mdl.joints));

mdl.calculateSensorJacobians();
dT10 = sens1.dTee_b;
dT02 = sens2.dTb_ee;
T10 = SE3.fastInverse(sens1.transform);
T02 = sens2.transform;

for i=1:numel(mdl.joints)
    dTee1_ee2(:,:,i) = dT10(:,:,i)*T02 + T10*dT02(:,:,i);
end

dT10_mul_T02 = reshape(permute(dT10,[2 1 3]),size(dT10,2),[])'*T02;
T10_mul_dT02 = (T10*reshape(permute(dT02,[2 1 3]),size(dT02,2),[]))';
dTee1_ee2 = permute(reshape((dT10_mul_T02 + T10_mul_dT02)',4,4,[]),[2 1 3]);




