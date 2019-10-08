

figNum = 1;

plotNames = {'p0','p1','p2','r0','r1','r2', ...
            'backFB','backAxial','backLateral', ...
            'rshldrPrism','rshldrElev','rshldrAbd','rshldrExtRot','relbowFlex','relbowSup', ...
            'lshldrPrism','lshldrElev','lshldrAbd','lshldrExtRot','lelbowFlex','lelbowSup', ...
            'rhipFlex','rhipAbd','rhipExtRot', ...
            'rkneeExtend','rkneeExtRot','rankleDorsi','ranklePron', ...
            'lhipFlex','lhipAbd','lhipExtRot', ...
            'lkneeExtend','lkneeExtRot','lankleDorsi','lanklePron'};
jumpNames = {'1','2','3','4','5','6','7','8','9','10','11','12'};
colorVec = [0,0.6,0; 0,0.9,0; 0.7,0.95,0; 1,0.9,0; 1,0.75,0; 1,0.55,0; 1,0.1,0.1; 0.95,0,0.6; ...
            0.4,0.1,0.7; 0,0.5,1; 0,0.9,0.95; 0,1,0.9];



% partsToLoad = {'02','03','04','05','06','07','08','09','10','11','12',...
%               '13','14','15','16','17','18','19','20','21','22'};
partsToLoad = {'10'};
          numParts = numel(partsToLoad);

for i_part = 1:numParts
    partNum = partsToLoad{i_part};

    % load JA data struct from certain participant
    mydir = pwd;
    idcs = strfind(mydir,'\');
    newdir = mydir(1:idcs(end)-1); 
    loadFilePath = [newdir '\results\JA_P' partNum '.mat'];
%     loadFilePath = [newdir '\results\JA_P' partNum '_2D.mat'];
    if(exist(loadFilePath,'file')~=2)
        disp(['JA_P' partNum ' file NOT FOUND']);
    else
        load(loadFilePath); % JA data struct, human data
        
        % (10-15) rshldrPrism, rshldrElev, rshldrAbd, rshldrExtRot, relbowFlex, relbowSup
        % (22-24) rhipFlex, rhipAbd, rhipExtRot, 
        % (25-28) rkneeExtend, rkneeExtRot, rankleDorsi, ranklePron
        
%         plotJoints = {[7,8],[10,15],[11,16],[12,17],[13,18],[20,27],[23,30],[25,32]};
%         plotJoints = {[11],[14],[22],[25],[27]};
        plotJoints = {[11,14,22,25,27]};
%         plotJA_multi(JA,plotJoints,1);
        
        
        
        for targSet = 1:2
            if(targSet==1)
                targCurr = JA.targ;
            else
                targCurr = JA.targAlign;
            end
        
            for i_targ = 2 %1:numel(targCurr)
                idx_LandFrame = JA.LandAlignJump(i_targ);
                currFig = figure(figNum); clf;
                currFig.Position = [50 50 800 600];
                
                TOFrame = JA.TOFrame(idx_LandFrame,i_targ);
                framesToPlot = (TOFrame-500):(TOFrame+500);
                
                
                for j1 = 1:numel(plotJoints)
                    clf; 
                    for j2 = 1:numel(plotJoints{j1})
                        subplot(numel(plotJoints{j1}),1,j2); hold on;
%                         sig = zeros(size(targCurr(i_targ).jump(1).data,1),numel(targCurr(i_targ).jump));
                        sig = zeros(numel(framesToPlot),numel(targCurr(i_targ).jump));

                        % get all jump data for specific joint
                        for i = 1:numel(targCurr(i_targ).jump)
                            sig(:,i) = targCurr(i_targ).jump(i).data(framesToPlot,plotJoints{j1}(j2));
                        end


                        % plot flipped/shifted data
                        for i = 1:size(sig,2)
                            if( sum(sig(:,i)) ~= 0 )
                                plot(rad2deg(sig(:,i)),'LineWidth',2,'Color',colorVec(i,:));
                            else %file wasn't available
                                plot(rad2deg(sig(:,i)),'b-.');
                            end
                        end
                        
                        if(targSet==2)
                            plot([JA.TOFrame(idx_LandFrame,i_targ), JA.TOFrame(idx_LandFrame,i_targ)]-(TOFrame-500),[-360,360],'k--');
                            plot([JA.LandFrame(idx_LandFrame,i_targ), JA.LandFrame(idx_LandFrame,i_targ)]-(TOFrame-500),[-360,360],'k');
                        end
%                         plot([JA.startFrame(idx_LandFrame,i_targ), JA.startFrame(idx_LandFrame,i_targ)],[-360,360],'k--');
%                         plot([JA.endFrame(idx_LandFrame,i_targ), JA.endFrame(idx_LandFrame,i_targ)],[-360,360],'k--');
                        
                        title(plotNames(plotJoints{j1}(j2)));
                        
                        xlim([0,1000]);
                        ylabel('Angle [deg]');
                        ylim([-170,170]);
                        grid on;
                        
                        if(j2 == numel(plotJoints{j1}))
                            legend(jumpNames(1:numel(targCurr(i_targ).jump)),'Location','northwest');
%                             xlabel('Frame Number');
                        end
                    end

        %             suptitle([eulTitle(eulOrLie) ', Target ' num2str(i_targ)]);
                    suptitle(['Part. ' num2str(JA.partNum) ': Target ' num2str(i_targ)]);
            %         ylabel('Angle [deg]');
                    xlabel('Frame Number');

        %             if(saveFigures)
        %                 plotName = plotNames(plotJoints{j1}(j2));
        %                 if(plotName
        %                 filename = ['P' num2str(JA.partNum) ': Target ' num2str(i_targ)]
        %                 saveas(gcf,filename);


%                     w = waitforbuttonpress;
                end

                figNum = figNum+1;
            end
        
        
        end
    end
end
        




