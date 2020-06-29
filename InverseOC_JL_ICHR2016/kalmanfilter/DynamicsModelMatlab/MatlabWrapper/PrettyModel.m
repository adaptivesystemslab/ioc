%This is just a file to visualize lower body model with a bunch of attached
%sensors

clear vis
clear mdl

video_path = 'D:\aslab\projects\vjoukov\General_FKEKF\DynamicsModelMatlab\MatlabWrapper\results\video';

%Load up XML model
mdl = rlCModel('Lower_Body_Rev_NoAnkle.xml');
mdl.forwardPosition();

%Create Visializer and place model in it
vis = rlVisualizer('vis',960,720);
vis.addModel(mdl);
vis.update();

%Change Visualizer Background to white 
vis.setBackgroundColour([1 1 1]);

%Make Video Writer and open it
vid = VideoWriter(video_path,'MPEG-4');
vid.open;

%Set Joint Angles Update Vis take screenshots and shove them into video
% KEEP VISUALIZER IN FOCUS OR IT WONT CAPTURE FRAMES
for i=1:1000
    
    %SET JOINT ANGLES HERE
    
    %THIS IS TO SET THEM IN ORDER
    mdl.position(:) = i/100;
    
    %YOU CAN ALSO ACCESS EACH JOINT SEPARATELY, THIS IS USEFUL IF YOU NEED
    %TO FIND THEM BY NAME (mdl.joints(2).name)
    mdl.joints(2).position = i/100;
    
    %MUST CALL THIS TO UPDATE MODEL
    mdl.forwardPosition();
    
    %UPDATE VIS AND CAPTURE FRAME
    vis.update();
    [vis_frame, b] = vis.getScreenshot();
    if b
       writeVideo(vid,vis_frame);
    end
end
%Dont forget to flush
vid.close();

%Here is how to make nice figure, copy pasted from other code so wont
%actually run.
fh2 = figure(1);
clf
fh2.Position = [900,100,720,480];
fh2.Color = [1 1 1];
title('EKF error during gimbal lock','FontSize',30,'FontWeight','bold');
ylabel('||I-R_{est}R_{gt}^T||_F','FontSize',25,'FontWeight','bold');
xlabel('Time (s)','FontSize',20,'FontWeight','bold')
set(gca,'fontsize',20,'fontweight','bold');

%THIS IS THE FANCY THING, USEFUL FOR VIDEO MAKING
ekf_line = animatedline('Color',[0.8500 0.3250 0.0980],'LineWidth',5);

%Open up 2 video writers, 1 for each figure
vid1 = VideoWriter('results\gimbal_lock_lg_ekf_err','MPEG-4');
vid1.open();

for k = 1:numel(t)

    %YOU ADD POINTS TO LINE AND CALL DRAWNOW 
    addpoints(ekf_line,t(k),ekf_rot_err(k));
    drawnow
    
    %THIS CAPTURES THE WHOLE FIGURE
    fh1_frame = getframe(fh1);
    
    %THEN PULL IMAGE OUT OF IT
    I1 = fh1_frame.cdata;
    
    writeVideo(vid1,I1);
end

%Close writers
vid1.close();





