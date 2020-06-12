%% This is a simple test to test rungeKuttaNystrom integration

%clear in correct order if cars exist
clear vis
clear mdl
clear all
%Create a visualizer
vis = Visualizer('vis',960,720);

%load up a model
mdl = rlCModel('planar.xml');
mdl.forwardPosition();


%Create 20 srnsors and add them to the model

for i=1:100
   sens =  SensorCore(['sensor' num2str(i)]);
   sens.addDecorator('accelerometer'); 
   sens.addDecorator('gyroscope'); 
   sens.addDecorator('position');
   sens.addDecorator('velocity');
   T1 = eye(4);
   T1(1:3,end) = (rand(3,1)-0.5);
   mdl.addSensor(sens,'frame23',T1);
end

mdl.forwardPosition();

%Add model to the visualizer
vis.addModel(mdl);
vis.update();
pause(10);

% [mdl.bodies.m] = deal(0.01);
% mdl.bodies(end).m = 0.1;

%Integrate with no input torques and see what happens
t = 0:0.01:100;
timer0 = zeros(size(t));
timer1 = zeros(size(t));
for i=1:numel(t)
   tic;
   mdl.rungeKuttaNystrom(0.01);
   mdl.forwardPosition();
   mdl.calculateSensorJacobians();
   %vis.update();
   timer0(i) = toc;
   
   %simulate friction on joints 
   mdl.torque = -mdl.velocity*0.2;
   
   %Pause dt
   pause(0.01-timer0(i));
   disp(i);
end

clear vis
clear mdl
clear sens1