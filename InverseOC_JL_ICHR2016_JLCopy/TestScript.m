%% Test Script
cd 'C:\Users\Rollen\Documents\General_FKEKF\DynamicsModelMatlab\MatlabWrapper'

vis = rlVisualizer('hello', 800, 600);
mdl = rlCModel('camera_arm.xml');
cam = Camera('MyCamera', mdl, 'framecamera', eye(4));

vis.addModel(mdl);
%vis.add3DModel('Box.3ds', [5 0 2 0 0 0]);
mdl.forwardPosition;

vis.addCamera(cam);
vis.update;

%%

% Decide on some aspect ratio
aspect = 16/9;

% Vertical Field Of View. The horizontal field of view is related by the
% aspect ratio.
vFOV = deg2rad(30);
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

% videoOutput = VideoWriter('C:\Users\Rollen\Videos\Features.mp4','MPEG-4');
% videoOutput.FrameRate = 24;
% videoOutput.open();
% figure;

I1 = [];
C1 = [];
V1 = [];
for t = 0:240
    mdl.joints(1).position = sin(2 * pi * t / 240) * (pi/16);
    mdl.joints(2).position = cos(2 * pi * t / 240) * (pi/16);
    mdl.forwardPosition

    % Set the pose of the camera and cache jacobians
    cam.updatePose;
    mdl.calculateSensorJacobians;
    
    % Update the visualizer
    vis.update;

    % Get an observation and cache the features.
    I0 = I1;
    C0 = C1;
    V0 = V1;
    I1 = cam.observation;
    [C1, V1, ~] = cam.getFeatures;
    
    if t == 0
    	continue;
    end
    
    % Print Jacobians for Cube Center
    cube_world = [4.5; 0; 2; 1];
    cube_pixel = cam.projectPoint(cube_world);
    cube_Jf = cam.jacobian_x0(cube_world);
    cube_Jq = cam.jacobian_q(cube_world);

%     
%     featureIndices = matchFeatures(C0, C1);
%     matchedV0 = V0(featureIndices(:,1),:);
%     matchedV1 = V1(featureIndices(:,2),:);
%     showMatchedFeatures(I0,I1,matchedV0,matchedV1);
%     imshow(I1);
%     hold on; 
%     plot(cube_pixel(1), cube_pixel(2), 'r.', 'MarkerSize', 20);
%     hold off;
%     drawnow
%     videoOutput.writeVideo(getframe);
end
% videoOutput.close();

%%
clear cam;
clear mdl;
clear vis;
clear all;

%%
clear DynamicsModelMatlab
exit