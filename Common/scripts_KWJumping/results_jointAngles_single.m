% Reads in previously-saved joint angle data and makes plots
clear; clc;
% Choose dataset and parameters

plotJointAngles = 1;


partNum = '02';
targetDist = '55';
i_set = 1;
i_jump = 6;

dataName = ['P' partNum '_target_' targetDist '_' num2str(i_set) ...
            '_' num2str(i_jump) '_clean-P' partNum '_template'];

mydir = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1); 
saveFilePath = [newdir '\results\RESULTS_JointAngles_P' partNum '_' targetDist '_' num2str(i_set) '_' num2str(i_jump) '.mat'];
load(saveFilePath);


if(plotJointAngles)
    %Plot joint angles
    % joints: (1-6) p0,p1,p2,r0,r1,r2, 
    % (7-9) backFB, backAxial, backLateral, 
    % (10-14) rshldrElev, rshldrAbd, rshldrExtRot, relbowFlex, relbowSup
    % (15-19) lshldrElev, lshldrAbd, lshldrExtRot, lelbowFlex, lelbowSup
    % (20-22) rhipFlex, rhipAbd, rhipExtRot, 
    % (23-26) rkneeExtend, rkneeExtRot, rankleDorsi, ranklePron
    % (27-29) lhipFlex, lhipAbd, lhipExtRot, 
    % (30-33) lkneeExtend, lkneeExtRot, lankleDorsi, lanklePron
    plotJoints(1).a = [7,8,9]; %back
    plotJoints(2).a = [10,15]; %shldr elevation
    plotJoints(3).a = [11,16]; %shldr abduction
    plotJoints(4).a = [12,17]; %shldr ext. rotation
    plotJoints(5).a = [13,18]; %elbow flexion
    plotJoints(6).a = [20,27]; %hip flexion
    % plotJoints(7).a = [21,28]; %hip abduction
    % plotJoints(8).a = [22,29]; %hip ext. rotation
    plotJoints(7).a = [23,30]; %knee flexion 
    % plotJoints(10).a = [24,31]; %knee ext. rotation
    plotJoints(8).a = [25,32]; %ankle dorsiflexion
    plotJoints(9).a = [26,33]; %ankle pronation
    plotNames = {'Back','Shoulder Elev.','Shoulder Abd.','Shoulder Ext. Rot.',...
        'Elbow Flex.','Hip Flex.','Knee Flex.','Ankle Dorsiflex.','Ankle Pronate'};

    dataLength = size(jointAngles_eul,1);
    figure(1);
    for i = 1:numel(plotJoints)
        clf; hold on;
        leg = [];
        for j = 1:numel(plotJoints(i).a)
            jj = plotJoints(i).a(j);
            plot(1:dataLength,rad2deg(jointAngles_eul(:,jj)),'LineWidth',2);
            leg = [leg; jointAngles_names(jj)];
        end
        plot([footTOFrame, footTOFrame],[-180,180],'g','LineWidth',2);
        plot([footLandFrame, footLandFrame],[-180,180],'r','LineWidth',2);
        plot([startFrame, startFrame],[-180,180],'k--');
        plot([endFrame, endFrame],[-180,180],'k--');
        title(plotNames(i));
        legend(cellstr(leg));
        xlabel('Frame Number');
        ylabel('Angle [deg]');
        ylim([-120,120]);
        grid on;
        w = waitforbuttonpress;
    end
end






