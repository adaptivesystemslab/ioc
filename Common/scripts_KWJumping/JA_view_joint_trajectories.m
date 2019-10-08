% load 

% partsToLoad = {'20'};
partsToLoad = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};
% partsToLoad = {'02','10','17','20'}; % "expert" jumpers, P02 is questionable
numParts = numel(partsToLoad);

plotSingleJump = 0;
plotSingleJoint = 1;
jointsToPlot = [4,5,6,9,10,11]; % back_jFB and right side arm and leg joints
jointNames = {'back-FB','rshldr','relbow','rhip','rknee','rankle'};
jointPlotYLim = [-2,2; -1,3; -1,3; -1,3; -3,1; -2,2];
numJoints = numel(jointsToPlot);


colorVec = [0.1,0.1,1; 0.8,0.5,1; 0,0.9,0; 0,0.9,0; 0.9,0.7,0; 0.9,0.1,0];
%        = [dark blue, light purple, green, green, yellow/orange, red]
lineVec = {'-','-','-','--','-','-'};
gradeTypes = {'B','SB','P','P*','SF','F'};
dt = 1/200;

figNum = 1;

if(plotSingleJump)
    flightTime = zeros(36,21);

    for i_part = 1:numParts
        partNum = partsToLoad{i_part};

        % load JA data struct from certain participant
        mydir = pwd;
        idcs = strfind(mydir,'\');
        newdir = mydir(1:idcs(end)-1); 
        loadFilePath = [newdir '\results\JA_P' partNum '_2D.mat'];
        if(exist(loadFilePath,'file')~=2)
            disp(['JA_P' partNum '_2D.mat NOT FOUND']);
        else
            load(loadFilePath); % JA data struct, human data

            partNum = JA.partNum;

            for i_targ = 1:3
                for i_jump = 1:12

                    jointTraj = JA.targ(i_targ).jump(i_jump).data(:,jointsToPlot);
                    time = [0:size(jointTraj,1)-1]*dt;

                    TOFrame = JA.TOFrame(i_jump,i_targ);
                    LandFrame = JA.LandFrame(i_jump,i_targ);

                    takeoffFrames = 200; % 1 second before takeoff frame ...
                    landFrames = 300; % ... to 1.5 seconds after takeoff frame (~ 1 second after landing)
                    framesToUse = (TOFrame - takeoffFrames):(TOFrame + landFrames); 

                    flightTime((12*(i_targ-1) + i_jump),i_part) = (LandFrame-TOFrame)*dt;

                    figure(figNum); clf; hold on; grid on;
                    for i_plot = 1:length(jointsToPlot)
                        plot(time, jointTraj(:,i_plot),'LineWidth',2);
                    end

                    plot([TOFrame*dt, TOFrame*dt],[-1.5,1.5],'k');
                    plot([LandFrame*dt, LandFrame*dt],[-1.5,1.5],'k');

                    xlim([(TOFrame - takeoffFrames)*dt,(TOFrame + landFrames)*dt]);
                    title(['P' partNum '-T' num2str(i_targ) '-J' num2str(i_jump)]); 
                    xlabel(['Flight time = ' num2str((LandFrame-TOFrame)*dt)]);

                    w = waitforbuttonpress();
                end
            end

        end  

    end
    
    figNum = figNum + 1;
end


if(plotSingleJoint)
    
    
    
    for i_part = 1:numParts
        partNum = partsToLoad{i_part};

        % load JA data struct from certain participant
        mydir = pwd;
        idcs = strfind(mydir,'\');
        newdir = mydir(1:idcs(end)-1); 
        loadFilePath = [newdir '\results\JA_P' partNum '_2D.mat'];
        if(exist(loadFilePath,'file')~=2)
            disp(['JA_P' partNum '_2D.mat NOT FOUND']);
        else
            load(loadFilePath); % JA data struct, human data
            partNum = JA.partNum;
            
            for i_targ = 1:3
                
                TOFrame = JA.TOFrame(JA.LandAlignJump(i_targ),i_targ);
                LandFrame = JA.LandFrame(JA.LandAlignJump(i_targ),i_targ);
                takeoffFrames = 200; % 1 second before takeoff frame ...
                landFrames = 300; % ... to 1.5 seconds after takeoff frame (~ 1 second after landing)
                framesToUse = (TOFrame - takeoffFrames):(TOFrame + landFrames);
                
                
                figure(figNum + i_targ - 1); clf;
                for i_plot = 1:numJoints
                    subplot(3,2,i_plot); hold on; grid on;
                    
                    for i_jump = 1:12
                        jointTraj = JA.targAlign(i_targ).jump(i_jump).data(:,jointsToPlot);
                        currGrade = char(JA.jumpGrades{12*(i_targ-1)+i_jump});
                        plotColorNum = find(strcmp(currGrade,gradeTypes)==1);
                        plot([1:size(jointTraj,1)]*dt, jointTraj(:,i_plot),'LineWidth',1,'Color',colorVec(plotColorNum,:));
                        
                    end
                    
                    plot([TOFrame*dt, TOFrame*dt],[-3,3],'k');
                    plot([LandFrame*dt, LandFrame*dt],[-3,3],'k');

                    xlim([(TOFrame - takeoffFrames)*dt,(TOFrame + landFrames)*dt]);
                    title(['P' partNum '-T' num2str(i_targ) '-Joint: ' jointNames{i_plot}]); 
                    xlabel('Time [sec]');
                    ylabel('Angle [rad]');
                    ylim(jointPlotYLim(i_plot,:));
                    
                end
                
            end
            
            w = waitforbuttonpress();
        end  

    end
    
    
end


close all;



