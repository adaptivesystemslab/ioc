% load 

% partsToLoad = {'03'};
partsToLoad = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};
% partsToLoad = {'02','10','17','20'}; % "expert" jumpers, P02 is questionable
numParts = numel(partsToLoad);

plotAllTarg = 0;
plotVelVec = 0;
plotIndivTarg_subplots = 1;
plotIndivTarg_fullPlots = 1;
plotVelVecOnly_subplots = 0;
plotVelVecOnly_fullPlots = 0;

plot_pos_vs_time_subplots = 0; % doesn't seem to show anything useful/differentiable
plot_vel_vs_time_subplots = 0; % doesn't seem to show anything useful/differentiable
plotAxes = [-1.3,0.8,-0.4,4.0,0.43,1.38];
pauseActive = 1;

savePlotsPNG = 1;
savePlotsFIG = 0;
saveFilePath = 'C:\Users\kgwester\Documents\ResearchWork\Thesis\Figures\CoM Trajectories\';
saveFilePathVelVec = 'C:\Users\kgwester\Documents\ResearchWork\Thesis\Figures\TOVel Plots\';

colorVec = [0.1,0.1,1; 0.8,0.5,1; 0,0.9,0; 0,0.9,0; 0.9,0.7,0; 0.9,0.1,0];
%        = [dark blue, light purple, green, green, yellow/orange, red]
lineVec = {'-','-','-','--','-','-'};
gradeTypes = {'B','SB','P','P*','SF','F'};
dt = 1/200;

figNum = 1;

for i_part = 1:numParts
    partNum = partsToLoad{i_part};
    
    % load JA data struct from certain participant
    mydir = pwd;
    idcs = strfind(mydir,'\');
    newdir = mydir(1:idcs(end)-1); 
    loadFilePath = [newdir '\results\JA_P' partNum '.mat'];
    if(exist(loadFilePath,'file')~=2)
        disp(['JA_P' partNum '.mat NOT FOUND']);
    else
        load(loadFilePath); % JA data struct, human data
%         jumpGrades = [JA.jumpGrades(1:3,:), JA.jumpGrades(4:6,:)];
        jumpGrades = JA.jumpGrades;
        
        if(savePlotsPNG || savePlotsFIG)
            figNum = 1;
        end

        %% Plot all CoM position and takeoff velocity vector on single plot
        if(plotAllTarg)
%                 grades = reshape(JA.jumpGrades',36,1);
            currFig = figure(figNum); clf; hold on; grid on;
            currFig.Position = [50 50 1000 800];
            xlabel('X'); ylabel('Y'); zlabel('Z');
            view([-30,30]);
            axis(plotAxes);

            plotOffset = 0;
            for currJump = 1:36
                trajEndFrame = find(squeeze(JA.CoMTraj(currJump,:,:))==[0,0,0],1,'first') - 30;
                 % "-30" just to get rid of appended zeros from shift alignment
                if(isempty(trajEndFrame))
                    traj = JA.CoMTraj(currJump,1:end-30,:);
                else
                    traj = JA.CoMTraj(currJump,1:trajEndFrame,:);
                end
                traj = squeeze(traj);

                i_targ = ceil(currJump/3);
                TOFrame = JA.TOFrame(JA.LandAlignJump(i_targ),i_targ);
                LandFrame = JA.LandFrame(JA.LandAlignJump(i_targ),i_targ);
                velPlotScale = 1/4.5; % scale vel vector to be visible in plot
                TOVelX = velPlotScale*(traj(TOFrame+1,1) - traj(TOFrame,1))/dt;
                TOVelZ = velPlotScale*(traj(TOFrame+1,3) - traj(TOFrame,3))/dt;

                currGrade = char(jumpGrades{currJump});
                plotColorNum = find(strcmp(currGrade,gradeTypes)==1);
                plot3(traj(:,1), traj(:,2)+plotOffset, traj(:,3),'Color',colorVec(plotColorNum,:),'LineStyle',lineVec{plotColorNum});

                if(mod(j,6)==0)
                    startLoc = JA.locationStart(currJump);
                    plot3([startLoc,startLoc],[plotOffset,plotOffset],[0,1],'k--');
                    targLoc = JA.locationLand(currJump,1);
                    plot3([targLoc(1),targLoc(1)],[plotOffset,plotOffset],[0,1],'k');
                    plotOffset = plotOffset + 0.1;
                end

            end


            figNum = figNum + 1;
        end

        %% Plot X vs. Z CoM position and takeoff velocity vector
        if(plotIndivTarg_subplots)
            currFig = figure(figNum); clf;
            currFig.Position = [50 50 1400 800];


            for i_targ = 1:3
                for i_set = 1:2
                    setNum = 2*(i_targ-1) + i_set;

                    ax(setNum) = subplot(3,2,setNum); hold on; grid on;
                    xlabel('X'); ylabel('Y'); zlabel('Z');
%                         view([-30,30]);
                    view([0,0]);
                    axis(plotAxes);

                    plotOffset = 0;
                    for j = 1:6
                        currJump = 12*(i_targ-1)+6*(i_set-1) + j;

                        trajEndFrame = find(squeeze(JA.CoMTraj(currJump,:,:))==[0,0,0],1,'first') - 30;
                         % "-30" just to get rid of any appended zeros from shift alignment
                        if(isempty(trajEndFrame))
                            traj = JA.CoMTraj(currJump,1:end-30,:);
                        else
                            traj = JA.CoMTraj(currJump,1:trajEndFrame,:);
                        end
                        traj = squeeze(traj);


                        TOFrame = JA.TOFrame(JA.LandAlignJump(i_targ),i_targ);
                        LandFrame = JA.LandFrame(JA.LandAlignJump(i_targ),i_targ);
                        velPlotScale = 1/4.5; % scale vel vector to be visible in plot
                        TOVelX = velPlotScale*(traj(TOFrame+1,1) - traj(TOFrame,1))/dt;
                        TOVelZ = velPlotScale*(traj(TOFrame+1,3) - traj(TOFrame,3))/dt;

                        startLoc = JA.locationStart(currJump);
                        velVecStart = [startLoc+0.25,0,0.5];

                        % Plot CoM traj and velocity vector based on
                        % jump grade
                        currGrade = char(JA.jumpGrades{currJump});
                        plotColorNum = find(strcmp(currGrade,gradeTypes)==1);

                        % CoM position traj, X vs. Z
                        plot3(traj(:,1), traj(:,2)+plotOffset, traj(:,3),'Color',colorVec(plotColorNum,:),'LineStyle',lineVec{plotColorNum});
                        % Takeoff CoM position
                        scatter3(traj(TOFrame,1), traj(TOFrame,2)+plotOffset, traj(TOFrame,3), 10, colorVec(plotColorNum,:),'filled');
                        scatter3(traj(TOFrame,1), traj(TOFrame,2)+plotOffset, traj(TOFrame,3), 10, 'k');
                        % Landing foot contact position
                        scatter3(traj(LandFrame,1), traj(LandFrame,2)+plotOffset, traj(LandFrame,3), 15, colorVec(plotColorNum,:),'filled');
                        scatter3(traj(LandFrame,1), traj(LandFrame,2)+plotOffset, traj(LandFrame,3), 15, 'k');
                        
                        % Takeoff velocity vector, origin at arbitrary
                        if(plotVelVec)
                            % aligned location
                            quiver3(velVecStart(1), velVecStart(2)+plotOffset, velVecStart(3), ...
                                TOVelX, 0, TOVelZ,'Color',colorVec(plotColorNum,:),'LineStyle',lineVec{plotColorNum}); 
                            % dot at end of velocity trajectory
                            scatter3(velVecStart(1)+TOVelX, velVecStart(2)+plotOffset, velVecStart(3)+TOVelZ, 4, colorVec(plotColorNum,:),'filled');
                        end

%                             startLoc = JA.locationStart(currJump);
                        plot3([startLoc,startLoc],[plotOffset,plotOffset],[0,2],'k--');
                        targLoc = JA.locationLand(currJump,1);
                        plot3([targLoc(1),targLoc(1)],[plotOffset,plotOffset],[0,2],'k');
                        plotOffset = plotOffset + 0.1;

                    end
                    title(['Part: ' partNum ', Target: ' num2str(i_targ) ', Set: ' num2str(i_set)]);
                    if(plotVelVec)
                        velVecPlotAxes(velVecStart,0);
                    end
                    
                end
            end

            hlink = linkprop(ax,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
            setappdata(gcf, 'StoreTheLink', Link);
%                 linkaxes;

            figNum = figNum + 1;

            if(savePlotsPNG)
                saveas(gcf,[saveFilePath 'Com_Traj_full_P' partNum '.png']);
            end
            if(savePlotsFIG)
                saveas(gcf,[saveFilePath 'Com_Traj_full_P' partNum '.fig']);
            end
        end


        %% Plot CoM position trajectories in individual full size plots 
        if(plotIndivTarg_fullPlots)
%                 figure(figNum); clf;

            lineThick = 1.5;

            for i_targ = 1:3
                for i_set = 1:2
                    currFig = figure(figNum); clf; hold on; grid on;
                    currFig.Position = [50 50 900 500];

                    xlabel('X'); ylabel('Y'); zlabel('Z');
                    view([0,0]);
                    axis(plotAxes);

                    plotOffset = 0;
                    for j = 1:6
                        currJump = 12*(i_targ-1)+6*(i_set-1) + j;

                        trajEndFrame = find(squeeze(JA.CoMTraj(currJump,:,:))==[0,0,0],1,'first') - 30;
                         % "-30" just to get rid of appended zeros from any shift alignment
                        if(isempty(trajEndFrame))
                            traj = JA.CoMTraj(currJump,1:end-30,:);
                        else
                            traj = JA.CoMTraj(currJump,1:trajEndFrame,:);
                        end
                        traj = squeeze(traj);

                        TOFrame = JA.TOFrame(JA.LandAlignJump(i_targ),i_targ);
                        LandFrame = JA.LandFrame(JA.LandAlignJump(i_targ),i_targ);
                        velPlotScale = 1/4.5; % scale vel vector to be visible in plot
                        TOVelX = velPlotScale*(traj(TOFrame+1,1) - traj(TOFrame,1))/dt;
                        TOVelZ = velPlotScale*(traj(TOFrame+1,3) - traj(TOFrame,3))/dt;

                        startLoc = JA.locationStart(currJump);
                        velVecStart = [startLoc+0.25,0,0.5];

                        % Plot CoM traj and velocity vector based on
                        % jump grade
                        currGrade = char(JA.jumpGrades{currJump});
                        plotColorNum = find(strcmp(currGrade,gradeTypes)==1);

                        % CoM position traj, X vs. Z
                        plot3(traj(:,1), traj(:,2)+plotOffset, traj(:,3),'Color',colorVec(plotColorNum,:),'LineStyle',lineVec{plotColorNum},'LineWidth',lineThick);
                        % Takeoff CoM position
                        scatter3(traj(TOFrame,1), traj(TOFrame,2)+plotOffset, traj(TOFrame,3), 20, colorVec(plotColorNum,:),'filled');
                        scatter3(traj(TOFrame,1), traj(TOFrame,2)+plotOffset, traj(TOFrame,3), 20, 'k');
                        % Landing foot contact position
                        scatter3(traj(LandFrame,1), traj(LandFrame,2)+plotOffset, traj(LandFrame,3), 30, colorVec(plotColorNum,:),'filled');
                        scatter3(traj(LandFrame,1), traj(LandFrame,2)+plotOffset, traj(LandFrame,3), 30, 'k');
                        
                        % Takeoff velocity vector, origin at arbitrary
                        if(plotVelVec)
                            % aligned location
                            quiver3(velVecStart(1), velVecStart(2)+plotOffset, velVecStart(3), ...
                                TOVelX, 0, TOVelZ,'Color',colorVec(plotColorNum,:),'LineStyle',lineVec{plotColorNum},'LineWidth',lineThick); 
                            % dot at end of velocity trajectory
                            scatter3(velVecStart(1)+TOVelX, velVecStart(2)+plotOffset, velVecStart(3)+TOVelZ, 10, colorVec(plotColorNum,:),'filled');
                        end
                        
%                             startLoc = JA.locationStart(currJump);
                        plot3([startLoc,startLoc],[plotOffset,plotOffset],[0,2],'k--');
                        targLoc = JA.locationLand(currJump,1);
                        plot3([targLoc(1),targLoc(1)],[plotOffset,plotOffset],[0,2],'k');
                        plotOffset = plotOffset + 0.1;

                    end
                    title(['Part: ' partNum ', Target: ' num2str(i_targ) ', Set: ' num2str(i_set)]);
                    if(plotVelVec)
                        velVecPlotAxes(velVecStart,1);
                    end

                    figNum = figNum + 1;

                    if(savePlotsPNG)
                        saveas(gcf,[saveFilePath 'Com_Traj_indiv_P' partNum '_T' num2str(i_targ) '_S' num2str(i_set) '.png']);
                    end
                    if(savePlotsFIG)
                        saveas(gcf,[saveFilePath 'Com_Traj_indiv_P' partNum '_T' num2str(i_targ) '_S' num2str(i_set) '.fig']);
                    end
                end
            end


        end
        
        
        
        
        
        
        %% Plot X vs. Z CoM position and takeoff velocity vector
        if(plotVelVecOnly_subplots)
            currFig = figure(figNum); clf;
            currFig.Position = [50 50 550 800];


            for i_targ = 1:3
                for i_set = 1:2
                    setNum = 2*(i_targ-1) + i_set;

                    ax(setNum) = subplot(3,2,setNum); hold on; grid on;
                    xlabel('$$\mathrm{\dot{X}}$$', 'Interpreter','latex'); 
                    ylabel('Y'); 
                    zlabel('$$\mathrm{\dot{Z}}$$', 'Interpreter','latex');
%                         view([-30,30]);
                    view([0,0]);
                    axis([0,2.5,0,2.5,0,2.5]);
                    

                    plotOffset = 0;
                    for j = 1:6
                        currJump = 12*(i_targ-1)+6*(i_set-1) + j;

                        trajEndFrame = find(squeeze(JA.CoMTraj(currJump,:,:))==[0,0,0],1,'first') - 30;
                         % "-30" just to get rid of any appended zeros from shift alignment
                        if(isempty(trajEndFrame))
                            traj = JA.CoMTraj(currJump,1:end-30,:);
                        else
                            traj = JA.CoMTraj(currJump,1:trajEndFrame,:);
                        end
                        traj = squeeze(traj);


                        TOFrame = JA.TOFrame(JA.LandAlignJump(i_targ),i_targ);
                        LandFrame = JA.LandFrame(JA.LandAlignJump(i_targ),i_targ);
                        velPlotScale = 1; % scale vel vector to be visible in plot
                        TOVelX = velPlotScale*(traj(TOFrame+1,1) - traj(TOFrame,1))/dt;
                        TOVelZ = velPlotScale*(traj(TOFrame+1,3) - traj(TOFrame,3))/dt;

%                         startLoc = JA.locationStart(currJump);
%                         velVecStart = [startLoc+0.25,0,0.5];
                        velVecStart = [0,0,0];

                        % Plot CoM traj and velocity vector based on
                        % jump grade
                        currGrade = char(JA.jumpGrades{currJump});
                        plotColorNum = find(strcmp(currGrade,gradeTypes)==1);

%                         % CoM position traj, X vs. Z
%                         plot3(traj(:,1), traj(:,2)+plotOffset, traj(:,3),'Color',colorVec(plotColorNum,:),'LineStyle',lineVec{plotColorNum});
%                         % Takeoff CoM position
%                         scatter3(traj(TOFrame,1), traj(TOFrame,2)+plotOffset, traj(TOFrame,3), 10, colorVec(plotColorNum,:),'filled');
%                         scatter3(traj(TOFrame,1), traj(TOFrame,2)+plotOffset, traj(TOFrame,3), 10, 'k');
%                         % Landing foot contact position
%                         scatter3(traj(LandFrame,1), traj(LandFrame,2)+plotOffset, traj(LandFrame,3), 15, colorVec(plotColorNum,:),'filled');
%                         scatter3(traj(LandFrame,1), traj(LandFrame,2)+plotOffset, traj(LandFrame,3), 15, 'k');
                        
                        % Takeoff velocity vector, origin at arbitrary
                        % aligned location
                        quiver3(velVecStart(1), velVecStart(2)+plotOffset, velVecStart(3), ...
                            TOVelX, 0, TOVelZ,'Color',colorVec(plotColorNum,:),'LineStyle',lineVec{plotColorNum}); 
                        % dot at end of velocity trajectory
                        scatter3(velVecStart(1)+TOVelX, velVecStart(2)+plotOffset, velVecStart(3)+TOVelZ, 4, colorVec(plotColorNum,:),'filled');

% %                             startLoc = JA.locationStart(currJump);
%                         plot3([startLoc,startLoc],[plotOffset,plotOffset],[0,2],'k--');
%                         targLoc = JA.locationLand(currJump,1);
%                         plot3([targLoc(1),targLoc(1)],[plotOffset,plotOffset],[0,2],'k');
%                         plotOffset = plotOffset + 0.1;

                    end
                    title(['Part: ' partNum ', Target: ' num2str(i_targ) ', Set: ' num2str(i_set)]);
%                     velVecPlotAxes(velVecStart,0);
                    axis([0,2.5,0,2.5,0,2.5]);
                    
                end
            end

%             hlink = linkprop(ax,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
%             setappdata(gcf, 'StoreTheLink', Link);
%                 linkaxes;
%             axis equal;
            
            figNum = figNum + 1;

            if(savePlotsPNG)
                saveas(gcf,[saveFilePathVelVec 'TOVel_full_P' partNum '.png']);
            end
            if(savePlotsFIG)
                saveas(gcf,[saveFilePathVelVec 'TOVel_full_P' partNum '.fig']);
            end
        end
        
        
        
        
        %% Plot CoM position trajectories in individual full size plots 
        if(plotVelVecOnly_fullPlots)
%                 figure(figNum); clf;

            lineThick = 1.5;

            for i_targ = 1:3
                for i_set = 1:2
                    currFig = figure(figNum); clf; hold on; grid on;
                    currFig.Position = [50 50 500 500];

                    xlabel('$$\mathrm{\dot{X}}$$', 'Interpreter','latex'); 
                    ylabel('Y'); 
                    zlabel('$$\mathrm{\dot{Z}}$$', 'Interpreter','latex');
                    view([0,0]);
                    axis([0,2.5,0,2.5,0,2.5]);

                    plotOffset = 0;
                    for j = 1:6
                        currJump = 12*(i_targ-1)+6*(i_set-1) + j;

                        trajEndFrame = find(squeeze(JA.CoMTraj(currJump,:,:))==[0,0,0],1,'first') - 30;
                         % "-30" just to get rid of appended zeros from any shift alignment
                        if(isempty(trajEndFrame))
                            traj = JA.CoMTraj(currJump,1:end-30,:);
                        else
                            traj = JA.CoMTraj(currJump,1:trajEndFrame,:);
                        end
                        traj = squeeze(traj);

                        TOFrame = JA.TOFrame(JA.LandAlignJump(i_targ),i_targ);
                        LandFrame = JA.LandFrame(JA.LandAlignJump(i_targ),i_targ);
                        velPlotScale = 1; % scale vel vector to be visible in plot
                        TOVelX = velPlotScale*(traj(TOFrame+1,1) - traj(TOFrame,1))/dt;
                        TOVelZ = velPlotScale*(traj(TOFrame+1,3) - traj(TOFrame,3))/dt;

                        startLoc = JA.locationStart(currJump);
%                         velVecStart = [startLoc+0.25,0,0.5];
                        velVecStart = [0,0,0];

                        % Plot CoM traj and velocity vector based on
                        % jump grade
                        currGrade = char(JA.jumpGrades{currJump});
                        plotColorNum = find(strcmp(currGrade,gradeTypes)==1);

%                         % CoM position traj, X vs. Z
%                         plot3(traj(:,1), traj(:,2)+plotOffset, traj(:,3),'Color',colorVec(plotColorNum,:),'LineStyle',lineVec{plotColorNum},'LineWidth',lineThick);
%                         % Takeoff CoM position
%                         scatter3(traj(TOFrame,1), traj(TOFrame,2)+plotOffset, traj(TOFrame,3), 20, colorVec(plotColorNum,:),'filled');
%                         scatter3(traj(TOFrame,1), traj(TOFrame,2)+plotOffset, traj(TOFrame,3), 20, 'k');
%                         % Landing foot contact position
%                         scatter3(traj(LandFrame,1), traj(LandFrame,2)+plotOffset, traj(LandFrame,3), 30, colorVec(plotColorNum,:),'filled');
%                         scatter3(traj(LandFrame,1), traj(LandFrame,2)+plotOffset, traj(LandFrame,3), 30, 'k');
                        
                        % Takeoff velocity vector, origin at arbitrary
                        % aligned location
                        quiver3(velVecStart(1), velVecStart(2)+plotOffset, velVecStart(3), ...
                            TOVelX, 0, TOVelZ,'Color',colorVec(plotColorNum,:),'LineStyle',lineVec{plotColorNum},'LineWidth',lineThick); 
                        % dot at end of velocity trajectory
                        scatter3(velVecStart(1)+TOVelX, velVecStart(2)+plotOffset, velVecStart(3)+TOVelZ, 10, colorVec(plotColorNum,:),'filled');
                        
% %                             startLoc = JA.locationStart(currJump);
%                         plot3([startLoc,startLoc],[plotOffset,plotOffset],[0,2],'k--');
%                         targLoc = JA.locationLand(currJump,1);
%                         plot3([targLoc(1),targLoc(1)],[plotOffset,plotOffset],[0,2],'k');
%                         plotOffset = plotOffset + 0.1;

                    end
                    title(['Part: ' partNum ', Target: ' num2str(i_targ) ', Set: ' num2str(i_set)]);
%                     velVecPlotAxes(velVecStart,1);
                    axis([0,2.5,0,2.5,0,2.5]);
                    
                    figNum = figNum + 1;

                    if(savePlotsPNG)
                        saveas(gcf,[saveFilePathVelVec 'TOVel_indiv_P' partNum '_T' num2str(i_targ) '_S' num2str(i_set) '.png']);
                    end
                    if(savePlotsFIG)
                        saveas(gcf,[saveFilePathVelVec 'TOVel_indiv_P' partNum '_T' num2str(i_targ) '_S' num2str(i_set) '.fig']);
                    end
                end
            end


        end
        
        
        
        
        %% Plot X vs. time and Z vs. time CoM position
        if(plot_pos_vs_time_subplots)
            currFig = figure(figNum); clf;
            currFig.Position = [50 50 1400 800];

            timeStart = 3.5; % don't plot first 3.5 seconds, CoM effectively doesn't move during calibration
            timeEnd = 7.5; % most people finish around this time
            for i_targ = 1:3
                for i_set = 1:2
                    setNum = 2*(i_targ-1) + i_set;

                    ax(setNum) = subplot(3,2,setNum); hold on; grid on;
                    xlabel('X'); ylabel('Y'); zlabel('Z');
%                         view([-30,30]);
                    view([0,0]);
                    axis([timeStart,timeEnd,-0.4,1.0,-1.7,1.7]);

                    plotOffset = 0;
                    for j = 1:6
                        currJump = 12*(i_targ-1)+6*(i_set-1) + j;

                        trajEndFrame = find(squeeze(JA.CoMTraj(currJump,:,:))==[0,0,0],1,'first') - 30;
                         % "-30" just to get rid of any appended zeros from shift alignment
                        if(isempty(trajEndFrame))
                            traj = JA.CoMTraj(currJump,1:end-30,:);
                        else
                            traj = JA.CoMTraj(currJump,1:trajEndFrame,:);
                        end
                        traj = squeeze(traj);


                        TOFrame = JA.TOFrame(JA.LandAlignJump(i_targ),i_targ);
                        LandFrame = JA.LandFrame(JA.LandAlignJump(i_targ),i_targ);
                        startLoc = JA.locationStart(currJump);
                        velVecStart = [startLoc+0.25,0,0.5];

                        % Plot CoM traj and velocity vector based on
                        % jump grade
                        currGrade = char(JA.jumpGrades{currJump});
                        plotColorNum = find(strcmp(currGrade,gradeTypes)==1);

                        % X pos vs. time
                        plot3((1:length(traj(:,1)))*dt, plotOffset*ones(length(traj(:,1)),1), traj(:,1),'Color',colorVec(plotColorNum,:),'LineStyle',lineVec{plotColorNum});
                        % Z pos vs. time
                        plot3((1:length(traj(:,3)))*dt, plotOffset*ones(length(traj(:,3)),1), traj(:,3),'Color',colorVec(plotColorNum,:),'LineStyle','-.');
                        % takeoff time
                        scatter3(TOFrame*dt, plotOffset, traj(TOFrame,1), 7, colorVec(plotColorNum,:),'filled');
                        scatter3(TOFrame*dt, plotOffset, traj(TOFrame,3), 7, colorVec(plotColorNum,:),'filled');
                        % landing time
                        scatter3(LandFrame*dt, plotOffset, traj(LandFrame,1), 9, colorVec(plotColorNum,:));
                        scatter3(LandFrame*dt, plotOffset, traj(LandFrame,3), 9, colorVec(plotColorNum,:));

%                             startLoc = JA.locationStart(currJump);
%                             plot3([startLoc,startLoc],[plotOffset,plotOffset],[0,2],'k--');
                        targLoc = JA.locationLand(currJump,1);
%                             plot3([targLoc(1),targLoc(1)],[plotOffset,plotOffset],[0,2],'k');
                        plotOffset = plotOffset + 0.1;

                    end
                    title(['Part: ' partNum ', Target: ' num2str(i_targ) ', Set: ' num2str(i_set)]);
                    xlabel('Time [sec]');


                end
            end

            hlink = linkprop(ax,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
            setappdata(gcf, 'StoreTheLink', Link);
%                 linkaxes;

            figNum = figNum + 1;
        end


        %% Plot X vs. time and Z vs. time CoM velocity
        if(plot_vel_vs_time_subplots)
            currFig = figure(figNum); clf;
            currFig.Position = [50 50 1400 800];

            timeStart = 3.5; % don't plot first 3.5 seconds, CoM effectively doesn't move during calibration
            timeEnd = 7.5; % most people finish around this time
            for i_targ = 1:3
                for i_set = 1:2
                    setNum = 2*(i_targ-1) + i_set;

                    ax(setNum) = subplot(3,2,setNum); hold on; grid on;
                    xlabel('X'); ylabel('Y'); zlabel('Z');
%                         view([-30,30]);
                    view([0,0]);
%                         axis([timeStart,timeEnd,-0.4,1.0,-1.7,1.7]);
                    xlim([timeStart, timeEnd]);

                    plotOffset = 0;
                    for j = 1:6
                        currJump = 12*(i_targ-1)+6*(i_set-1) + j;

                        trajEndFrame = find(squeeze(JA.CoMTraj(currJump,:,:))==[0,0,0],1,'first') - 30;
                         % "-30" just to get rid of any appended zeros from shift alignment
                        if(isempty(trajEndFrame))
                            traj = JA.CoMTraj(currJump,1:end-30,:);
                        else
                            traj = JA.CoMTraj(currJump,1:trajEndFrame,:);
                        end
                        traj = squeeze(traj);

                        % Compute velocity trajectory
                        trajVel = (traj(2:end,:)-traj(1:end-1,:))/dt;

                        TOFrame = JA.TOFrame(JA.LandAlignJump(i_targ),i_targ);
                        LandFrame = JA.LandFrame(JA.LandAlignJump(i_targ),i_targ);
                        startLoc = JA.locationStart(currJump);
                        velVecStart = [startLoc+0.25,0,0.5];

                        % Plot CoM traj and velocity vector based on
                        % jump grade
                        currGrade = char(JA.jumpGrades{currJump});
                        plotColorNum = find(strcmp(currGrade,gradeTypes)==1);

                        % X pos vs. time
                        plot3((1:length(trajVel(:,1)))*dt, plotOffset*ones(length(trajVel(:,1)),1), trajVel(:,1),'Color',colorVec(plotColorNum,:),'LineStyle',lineVec{plotColorNum});
                        % Z pos vs. time
                        plot3((1:length(trajVel(:,3)))*dt, plotOffset*ones(length(trajVel(:,3)),1), trajVel(:,3),'Color',colorVec(plotColorNum,:),'LineStyle','-.');
                        % takeoff time
                        scatter3(TOFrame*dt, plotOffset, trajVel(TOFrame,1), 7, colorVec(plotColorNum,:),'filled');
                        scatter3(TOFrame*dt, plotOffset, trajVel(TOFrame,3), 7, colorVec(plotColorNum,:),'filled');
                        % landing time
                        scatter3(LandFrame*dt, plotOffset, trajVel(LandFrame,1), 9, colorVec(plotColorNum,:));
                        scatter3(LandFrame*dt, plotOffset, trajVel(LandFrame,3), 9, colorVec(plotColorNum,:));

%                             startLoc = JA.locationStart(currJump);
%                             plot3([startLoc,startLoc],[plotOffset,plotOffset],[0,2],'k--');
                        targLoc = JA.locationLand(currJump,1);
%                             plot3([targLoc(1),targLoc(1)],[plotOffset,plotOffset],[0,2],'k');
                        plotOffset = plotOffset + 0.1;

                    end
                    title(['Part: ' partNum ', Target: ' num2str(i_targ) ', Set: ' num2str(i_set)]);
                    xlabel('Time [sec]');


                end
            end

            hlink = linkprop(ax,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
            setappdata(gcf, 'StoreTheLink', Link);
%                 linkaxes;

            figNum = figNum + 1;
        end



            
            
            
    end
%     w = waitforbuttonpress();   
    
end
    




