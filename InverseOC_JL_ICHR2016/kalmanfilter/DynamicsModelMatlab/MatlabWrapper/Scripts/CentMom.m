%% Test centroidal momentum stuff

%Load up a model visualize and compute centroidal momentum

%Add main folder
addpath('..\');
%Create Model and visualize
mdl = rlCModel('simple_3dof_model.xml');
mdl.forwardPosition();

vis = rlVisualizer('vis',640,480);
vis.addModel(mdl);
vis.update();

%% Compute Momentum in base frame

mdl.velocity = rand(size(mdl.velocity));
mdl.forwardVelocity();

X_star = zeros(6,6,numel(mdl.bodies));
Hg = zeros(6,1);
Hg_dot = zeros(6,1);
for i=1:numel(mdl.bodies)
    b = mdl.bodies(i);
    X = b.x;
    X_star(:,:,i) = X;
    X_star(1:3,4:6,i) = X(4:6,1:3);
    X_star(4:6,1:3,i) = 0;
    
    Hg = Hg + X_star(:,:,i)*b.i*b.v;
    
    % Velocity in ground frame
    Vg = b.x*b.v;
    % Cross product matrix of Vg
    V_star = [skew(Vg(1:3)) skew(Vg(4:6)); zeros(3,3) skew(Vg(1:3))];
    Hg_dot = Hg_dot + V_star*X_star(:,:,i)*b.i*b.v + X_star(:,:,i)*b.i*b.a;
end


