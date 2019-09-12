grades = reshape(JA.jumpGrades',36,1);
figure(1); clf; hold on; grid on;
for j = 1:36
    trajEndFrame = find(squeeze(JA.CoMTraj(j,:,:))==[0,0,0],1,'first') - 30;
     % "-30" just to get rid of appended zeros from shift alignment
    if(isempty(trajEndFrame))
        traj = JA.CoMTraj(j,1:end-30,:);
    else
        traj = JA.CoMTraj(j,1:trajEndFrame,:);
    end
    traj = squeeze(traj);
    
    if(contains(grades{j},'*'))
        plot3(traj(:,1), traj(:,2), traj(:,3),'r--');
    elseif(contains(grades{j},'P'))
        plot3(traj(:,1), traj(:,2), traj(:,3),'r');
    elseif(contains(grades{j},'S'))
        plot3(traj(:,1), traj(:,2), traj(:,3),'b');
    else
        plot3(traj(:,1), traj(:,2), traj(:,3),'b--');
    end
    axis([-1.5,1,-1,1,0.4,1.3]);
end

