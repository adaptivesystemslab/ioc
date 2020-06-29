%% Test Script
cd 'C:\Users\Rollen\Documents\General_FKEKF\DynamicsModelMatlab\MatlabWrapper'

vis = rlVisualizer('hello', 800, 600);
mdl = rlCModel('camera_arm.xml');
cam = Camera('MyCamera', mdl, 'framecamera', eye(4));
mdl.forwardPosition;

vis.add3DModel('Box.3ds', [5; 0; 2; 0; 0; 0;]);
vis.addModel(mdl);
vis.addCamera(cam);
vis.update;

%%

% Decide on some aspect ratio
aspect = 1;

% Vertical Field Of View. The horizontal field of view is related by the
% aspect ratio.
vFOV = deg2rad(90);
hFOV = 2 * atan(tan(vFOV / 2) * aspect);

% Width and Height of the Image.
width = 640;
height = width / aspect;

% Focal lengths (in pixel coordinates) are related to the image size and
% the field of view in the direction of concern. This is what is found in
% the OpenCV Calibration Matrix.
fX = (1 / tan(hFOV / 2)) * width / 2;
fY = (1 / tan(vFOV / 2)) * height / 2;

% Calibration Camera :)
cam.calibration = [fX,0,width/2;0,fY,height/2;0,0,1];

%%
data = [];
[Measurements, PActual] = testMeasurements(vis, mdl, cam);
sum(isnan(Measurements))

%%
mdl.joints(1).position = sin(2 * pi * 0 / 240) * (pi/32);
mdl.joints(2).position = cos(2 * pi * 0 / 53) * (pi/16);
mdl.joints(3).position = tanh((0 - 120)/120) * pi/16;
mdl.forwardPosition;

ekf = Camera_EKF_Q(mdl, cam);
MEstimate = [];
PEstimate = zeros(240, numel(mdl.joints));
for t = 1:240
    ekf.run_iteration(0.05, Measurements(t,:)');
    est = ekf.makeMeasure(ekf.state)';
    MEstimate = [MEstimate; est];
    PEstimate(t,:) = mdl.position;
%     cam.updatePose;
% 
%     vis.update;
% 
%     I = cam.observation;
%     imshow(I);
%     hold on;
%     for fi = 1:(length(est)/2)
%         plot(est(2*fi-1), est(2*fi), 'b+', 'MarkerSize', 15);
%     end
%     for fi = 1:(length(est)/2)
%         plot(Measurements(t,2*fi-1), Measurements(t,2*fi), 'r.', 'MarkerSize', 15);
%     end
%     hold off;
%     drawnow;
%     pause(0.05);
end

%% Error Calculation
dm = (MEstimate - Measurements) .^ 2;
error = zeros(size(dm,1), (size(dm, 2)/2));
for i=1:(size(dm, 2)/2)
    for t = 1:size(dm,1)
        error(t,i) = sqrt(sum(dm(t,(2*i-1):(2*i))));
    end
end
R = max(error)

%%
clear cam;
clear mdl;
clear vis;