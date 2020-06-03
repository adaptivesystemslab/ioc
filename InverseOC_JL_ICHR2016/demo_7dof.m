%7dof demo

%Clear in correct order
clear vis;
clear model;

%Load the model
model = rlCModel('7dof_squat.xml');
%Run forward position
model.forwardPosition();

%Create visualizer window
vis = rlVisualizer('vis',640,480);
%Add model to visualizer
vis.addModel(model);
%Update scene
vis.update();
