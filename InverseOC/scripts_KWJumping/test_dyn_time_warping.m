% function plotJA_multi(JA,plotJoints,figNum)

partNum = '03';
plotJoints = {[7],[10,15],[13,18],[20,27],[23,30],[25,32]}; %backFB, and all arm/leg joints in sagittal plane
% plotJoints = {[7],[10],[13],[20],[23],[25]}; %right sides only

figNum = 1;

% load JA struct
mydir = pwd;
idcs = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1); 
loadFilePath = [newdir '\results\JA_P' partNum '.mat'];
if(exist(loadFilePath,'file')~=2)
    disp(['JA_P' partNum '.mat NOT FOUND']);
else
    load(loadFilePath);
    targCurr = JA.targAlign_eul;

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
    
    shift_rec = zeros(11,numel(targCurr(i_targ).jump),numel(targCurr)); %"11" is number of joints to look at
    
    
    for i_targ = 1:numel(targCurr)
        idx_LandFrame = JA.LandAlignJump(i_targ);
        jointNum = 1;

        for j1 = 1:numel(plotJoints)
            figure(figNum); clf; 
            for j2 = 1:numel(plotJoints{j1})
                subplot(numel(plotJoints{j1}),2,2*j2-1); hold on;
                sig = zeros(size(targCurr(i_targ).jump(1).data,1),numel(targCurr(i_targ).jump));
                
                % get all jump data for specific joint
                for i = 1:numel(targCurr(i_targ).jump)
                    sig(:,i) = targCurr(i_targ).jump(i).data(:,plotJoints{j1}(j2));
                end
                
                % plot data
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
                legend(jumpNames(1:numel(targCurr(i_targ).jump)));
                title(plotNames(plotJoints{j1}(j2)));
    %             ylim([ymin,ymax]);
                ylim([-150,150]);
                grid on;
                
                
                

%%%%%%%%%%%%%% perform shifting on each signal
                subplot(numel(plotJoints{j1}),2,2*j2); hold on;
                sig_dtw = sig;
                
                %NOTE: only perform DTW on frames 100 (0.5 sec.) before TO
                % and 100 after landing
                dtw_st = JA.TOFrame(idx_LandFrame,i_targ) - 100;
                dtw_end = JA.LandFrame(idx_LandFrame,i_targ) + 100;
                maxsamp = 50; %max samples shifted by in DTW algorithm
                
                for i = 1:size(sig_dtw,2)-1 %use last jump as DTW reference
                    tmp = sig_dtw(dtw_st:dtw_end,i);
                    [dist,ix,iy] = dtw(sig_dtw(dtw_st:dtw_end,1),tmp,maxsamp);
                    
                    sig_dtw(:,i) = tmp(iy); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OR SOMETHING
                    
                    
                end
                
                
                
                % plot data
                for i = 1:size(sig_dtw,2)
                    if( sum(sig_dtw(:,i)) ~= 0 )
                        plot(rad2deg(sig_dtw(:,i)),'LineWidth',2,'Color',colorVec(i,:));
                    else %file wasn't available
                        plot(rad2deg(sig_dtw(:,i)),'b-.');
                    end
                end
                
                plot([JA.TOFrame(idx_LandFrame,i_targ), JA.TOFrame(idx_LandFrame,i_targ)],[-360,360],'g--');
                plot([JA.LandFrame(idx_LandFrame,i_targ), JA.LandFrame(idx_LandFrame,i_targ)],[-360,360],'r--');
                plot([JA.startFrame(idx_LandFrame,i_targ), JA.startFrame(idx_LandFrame,i_targ)],[-360,360],'k--');
                plot([JA.endFrame(idx_LandFrame,i_targ), JA.endFrame(idx_LandFrame,i_targ)],[-360,360],'k--');
                legend(jumpNames(1:numel(targCurr(i_targ).jump)));
                title(plotNames(plotJoints{j1}(j2)));
    %             ylim([ymin,ymax]);
                ylim([-150,150]);
                grid on;
                
                
                jointNum = jointNum + 1;
            end
            
            suptitle([num2str(JA.partNum) ' Target ' num2str(i_targ)]);
    %         ylabel('Angle [deg]');
            xlabel('Frame');
%             w = waitforbuttonpress;
        end

        figNum = figNum+1;
    end
end
        

