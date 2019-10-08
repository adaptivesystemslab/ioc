% load 


partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};
partsToRun = {'11'};




for i_part = 1:numel(partsToRun)
    partNum = partsToRun{i_part};
    
    % load JA data struct from certain participant
    mydir = pwd;
    idcs = strfind(mydir,'\');
    newdir = mydir(1:idcs(end)-1); 
    loadFilePath = [newdir '\results\JA_P' partNum '.mat'];
    if(exist(loadFilePath,'file')~=2)
        disp(['JA_P' partNum '.mat NOT FOUND']);
    else
        load(loadFilePath); % JA data struct, human data
        jumpGrades = [JA.jumpGrades(1:3,:), JA.jumpGrades(4:6,:)];
        
        for i_targ = 1:numel(JA.targAlign)
            figNum = 1;
            
            if(plotAllTarg)
                grades = reshape(JA.jumpGrades',36,1);
                currFig = figure(figNum); clf; hold on; grid on;
                currFig.Position = [50 50 1000 800];
                xlabel('X'); ylabel('Y'); zlabel('Z');
                view([-30,30]);
                axis([-1.5,1,-0.4,1.4,0.4,1.3]);
                
                plotOffset = 0;
                for j = 1:36
                    trajEndFrame = find(squeeze(JA.CoMTraj(j,:,:))==[0,0,0],1,'first') - 30;
                     % "-30" just to get rid of appended zeros from shift alignment
                    if(isempty(trajEndFrame))
                        traj = JA.CoMTraj(j,1:end-30,:);
                    else
                        traj = JA.CoMTraj(j,1:trajEndFrame,:);
                    end
                    traj = squeeze(traj);

                    currGrade = char(grades{j});
                    switch currGrade
                        case 'B'
                            plot3(traj(:,1), traj(:,2)+plotOffset, traj(:,3),'Color',colorVec(1,:));
                        case 'SB'
                            plot3(traj(:,1), traj(:,2)+plotOffset, traj(:,3),'Color',colorVec(2,:));
                        case 'P'
                            plot3(traj(:,1), traj(:,2)+plotOffset, traj(:,3),'Color',colorVec(3,:));
                        case 'P*'
                            plot3(traj(:,1), traj(:,2)+plotOffset, traj(:,3),'Color',colorVec(3,:),'LineStyle','--');
                        case 'SF'
                            plot3(traj(:,1), traj(:,2)+plotOffset, traj(:,3),'Color',colorVec(4,:));
                        case 'F'
                            plot3(traj(:,1), traj(:,2)+plotOffset, traj(:,3),'Color',colorVec(5,:));
                    end
                    
                    if(mod(j,6)==0)
                        startLoc = JA.locationStart(j);
                        plot3([startLoc,startLoc],[plotOffset,plotOffset],[0,1],'k--');
                        targLoc = JA.locationLand(j,:);
                        plot3([targLoc(1),targLoc(1)],[plotOffset,plotOffset],[0,1],'k');
                        plotOffset = plotOffset + 0.2;
                    end
                    
                end
                
                
                figNum = figNum + 1;
            end
            
            if(plotIndivTarg_subplots)
                currFig = figure(figNum); clf;
                currFig.Position = [50 50 1000 800];
                
                for i_targ = 1:3
                    for i_set = 1:2
                        setNum = 2*(i_targ-1) + i_set;
                        
                        grades = JA.jumpGrades{setNum,:};
                        subplot(3,2,setNum); hold on; grid on;
                        xlabel('X'); ylabel('Y'); zlabel('Z');
                        view([-30,30]);
                        axis([-1.5,1,-0.4,1.4,0.4,1.3]);

                        for j = 1:6
                            trajEndFrame = find(squeeze(JA.CoMTraj(j,:,:))==[0,0,0],1,'first') - 30;
                             % "-30" just to get rid of appended zeros from shift alignment
                            if(isempty(trajEndFrame))
                                traj = JA.CoMTraj(j,1:end-30,:);
                            else
                                traj = JA.CoMTraj(j,1:trajEndFrame,:);
                            end
                            traj = squeeze(traj);

                            currGrade = char(JA.jumpGrades{setNum,j});
                            switch currGrade
                                case 'B'
                                    plot3(traj(:,1), traj(:,2), traj(:,3),'Color',colorVec(1,:));
                                case 'SB'
                                    plot3(traj(:,1), traj(:,2), traj(:,3),'Color',colorVec(2,:));
                                case 'P'
                                    plot3(traj(:,1), traj(:,2), traj(:,3),'Color',colorVec(3,:));
                                case 'P*'
                                    plot3(traj(:,1), traj(:,2), traj(:,3),'Color',colorVec(3,:),'LineStyle','--');
                                case 'SF'
                                    plot3(traj(:,1), traj(:,2), traj(:,3),'Color',colorVec(4,:));
                                case 'F'
                                    plot3(traj(:,1), traj(:,2), traj(:,3),'Color',colorVec(5,:));
                            end

                            startLoc = JA.locationStart(6*setNum);
                            plot3([startLoc,startLoc],[plotOffset,plotOffset],[0,2],'k--');
                            targLoc = JA.locationLand(6*setNum,:);
                            plot3([targLoc(1),targLoc(1)],[plotOffset,plotOffset],[0,2],'k');
                            plotOffset = plotOffset + 0.2;

                        end
                        title(['Part: ' partNum ', Target: ' num2str(i_targ) ', Set: ' num2str(i_set)]);
                        
                        figNum = figNum + 1;
                    end
                end
                
                
            end
            
            
            if(plotIndivTarg_fullPlots)
%                 figure(figNum); clf;
                
                for i_targ = 1:3
                    for i_set = 1:2
                        currFig = figure(figNum); clf; hold on; grid on;
                        currFig.Position = [50 50 1000 800];
                        
                        setNum = 2*(i_targ-1) + i_set;
                        
                        grades = JA.jumpGrades{setNum,:};
                        xlabel('X'); ylabel('Y'); %zlabel('Z');
%                         view([-30,30]);
%                         axis([-1.5,1,-0.4,1.4,0.4,1.3]);
                        axis([-1.5,1,0.6,1.3]);
                        
                        for j = 1:6
                            trajEndFrame = find(squeeze(JA.CoMTraj(j,:,:))==[0,0,0],1,'first') - 30;
                             % "-30" just to get rid of appended zeros from any shift alignment
                            if(isempty(trajEndFrame))
                                traj = JA.CoMTraj(6*(setNum-1)+j,1:end-30,:);
                            else
                                traj = JA.CoMTraj(6*(setNum-1)+j,1:trajEndFrame,:);
                            end
                            traj = squeeze(traj);

                            currGrade = char(JA.jumpGrades{setNum,j});
                            switch currGrade
                                case 'B'
                                    plot(traj(:,1), traj(:,3),'Color',colorVec(1,:),'LineWidth',2);
                                case 'SB'
                                    plot(traj(:,1), traj(:,3),'Color',colorVec(2,:),'LineWidth',2);
                                case 'P'
                                    plot(traj(:,1), traj(:,3),'Color',colorVec(3,:),'LineWidth',2);
                                case 'P*'
                                    plot(traj(:,1), traj(:,3),'Color',colorVec(3,:),'LineWidth',2,'LineStyle','--');
                                case 'SF'
                                    plot(traj(:,1), traj(:,3),'Color',colorVec(4,:),'LineWidth',2);
                                case 'F'
                                    plot(traj(:,1), traj(:,3),'Color',colorVec(5,:),'LineWidth',2);
                            end

                            startLoc = JA.locationStart(6*setNum);
                            plot([startLoc,startLoc],[0,2],'k--');
                            targLoc = JA.locationLand(6*setNum,:);
                            plot([targLoc(1),targLoc(1)],[0,2],'k');
%                             plotOffset = plotOffset + 0.2;

                        end
                        title(['Part: ' partNum ', Target: ' num2str(i_targ) ', Set: ' num2str(i_set)]);
                        
                        figNum = figNum + 1;
                    end
                end
                
                
            end
            
            
%             w = waitforbuttonpress();
            
        end
    end
end



