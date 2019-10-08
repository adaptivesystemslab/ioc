function plotJA_multi_compare_shiftJAs(JA,plotJoints,figNum)
% "JA" is joint angle data struct

%Plot joints legend
% (1-6) p0,p1,p2,r0,r1,r2, 
% (7-9) backFB, backAxial, backLateral, 
% (10-14) rshldrElev, rshldrAbd, rshldrExtRot, relbowFlex, relbowSup
% (15-19) lshldrElev, lshldrAbd, lshldrExtRot, lelbowFlex, lelbowSup
% (20-22) rhipFlex, rhipAbd, rhipExtRot, 
% (23-26) rkneeExtend, rkneeExtRot, rankleDorsi, ranklePron
% (27-29) lhipFlex, lhipAbd, lhipExtRot, 
% (30-33) lkneeExtend, lkneeExtRot, lankleDorsi, lanklePron

% Example "plotJoints" input format:
% plotJoints = {[7,8,9],[10,15],[11,16],[12,17],[13,18],[20,27],[23,30],[25,32],[26,33]};

% "figNum" = number of first figure

plotNames = {'p0','p1','p2','r0','r1','r2', ...
            'backFB','backAxial','backLateral', ...
            'rshldrElev','rshldrAbd','rshldrExtRot','relbowFlex','relbowSup', ...
            'lshldrElev','lshldrAbd','lshldrExtRot','lelbowFlex','lelbowSup', ...
            'rhipFlex','rhipAbd','rhipExtRot', ...
            'rkneeExtend','rkneeExtRot','rankleDorsi','ranklePron', ...
            'lhipFlex','lhipAbd','lhipExtRot', ...
            'lkneeExtend','lkneeExtRot','lankleDorsi','lanklePron'};
jumpNames = {'1','2','3','4','5','6','7','8','9','10','11','12'};
colorVec = [0,0.6,0; 0,0.9,0; 0.7,0.95,0; 1,0.9,0; 1,0.75,0; 1,0.55,0; 1,0.1,0.1; 0.95,0,0.6; ...
            0.4,0.1,0.7; 0,0.5,1; 0,0.9,0.95; 0,1,0.9];

eulTitle = {'Eul','Lie'};
for eulOrLie = 1%:2 %add in when want to see lie joint angles too
    if(eulOrLie==1)
        targCurr = JA.targAlign_eul;
        targCurr_fine = JA.targAlign_eul_fine;
    else
        targCurr = JA.targAlign_lie;
    end
    
    
    for i_targ = 1:numel(targCurr)
        idx_LandFrame = JA.LandAlignJump(i_targ);
        currFig = figure(figNum); clf;
        currFig.Position = [100 100 600 800];
        
        for j1 = 1:numel(plotJoints)
            clf; 
            for j2 = 1:numel(plotJoints{j1})
                subplot(numel(plotJoints{j1}),2,2*j2-1); hold on;
                sig = zeros(size(targCurr(i_targ).jump(1).data,1),numel(targCurr(i_targ).jump));

                % get all jump data for specific joint
                for i = 1:numel(targCurr(i_targ).jump)
                    sig(:,i) = targCurr(i_targ).jump(i).data(:,plotJoints{j1}(j2));
                end


                % plot flipped/shifted data
                for i = 1:size(sig,2)
                    if( sig(:,i) ~= zeros(size(sig,1),1) )
                        plot(rad2deg(sig(:,i)),'LineWidth',2,'Color',colorVec(i,:));
                    else %file wasn't available
                        plot(rad2deg(sig(:,i)),'b-.');
                    end
                end

                plot([JA.TOFrame(idx_LandFrame,i_targ), JA.TOFrame(idx_LandFrame,i_targ)],[-360,360],'g--');
                plot([JA.LandFrame(idx_LandFrame,i_targ), JA.LandFrame(idx_LandFrame,i_targ)],[-360,360],'r--');
                plot([JA.startFrame(idx_LandFrame,i_targ), JA.startFrame(idx_LandFrame,i_targ)],[-360,360],'k--');
                plot([JA.endFrame(idx_LandFrame,i_targ), JA.endFrame(idx_LandFrame,i_targ)],[-360,360],'k--');
                legend(jumpNames(1:numel(targCurr(i_targ).jump)),'Location','northwest');
                title(plotNames(plotJoints{j1}(j2)));
                ylim([-170,170]);
                grid on;
%             end
            
%           
            
            
            
            %% second set of subplots, 
                subplot(numel(plotJoints{j1}),2,2*j2); hold on;
                sig = zeros(size(targCurr_fine(i_targ).jump(1).data,1),numel(targCurr_fine(i_targ).jump));

                % get all jump data for specific joint
                for i = 1:numel(targCurr_fine(i_targ).jump)
                    sig(:,i) = targCurr_fine(i_targ).jump(i).data(:,plotJoints{j1}(j2));
                end


                % plot flipped/shifted data
                for i = 1:size(sig,2)
                    if( sum(sig(:,i)) ~= 0 )
                        plot(rad2deg(sig(:,i)),'LineWidth',2,'Color',colorVec(i,:));
                    else %file wasn't available
                        plot(rad2deg(sig(:,i)),'b-.');
                    end
                end

                plot([JA.TOFrame(idx_LandFrame,i_targ), JA.TOFrame(idx_LandFrame,i_targ)],[-360,360],'g--');
                plot([JA.LandFrame(idx_LandFrame,i_targ), JA.LandFrame(idx_LandFrame,i_targ)],[-360,360],'r--');
                plot([JA.startFrame(idx_LandFrame,i_targ), JA.startFrame(idx_LandFrame,i_targ)],[-360,360],'k--');
                plot([JA.endFrame(idx_LandFrame,i_targ), JA.endFrame(idx_LandFrame,i_targ)],[-360,360],'k--');
                legend(jumpNames(1:numel(targCurr(i_targ).jump)),'Location','northwest');
                title(['fine align']);
                ylim([-170,170]);
                grid on;
            end
            linkaxes;
%             suptitle([eulTitle(eulOrLie) ', Target ' num2str(i_targ)]);
            suptitle(['Part. ' num2str(JA.partNum) ': Target ' num2str(i_targ)]);
    %         ylabel('Angle [deg]');
            
            
            
            w = waitforbuttonpress;
        end

        figNum = figNum+1;
    end
end
        

