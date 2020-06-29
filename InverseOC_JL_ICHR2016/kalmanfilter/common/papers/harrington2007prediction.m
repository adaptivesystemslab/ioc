function [global2RHC, global2LHC, mid2RHC, mid2LHC] = harrington2007prediction(trcData, cropInd1, cropInd2)

    %Rotation calculate based on the mean of mid->front and right for y
    %Here we get the world to base transform to match Mocap Data
    front = mean((trcData.ASIS_R(cropInd1:cropInd2,:)+trcData.ASIS_L(cropInd1:cropInd2,:))/2);
    back = mean((trcData.PSIS_R(cropInd1:cropInd2,:)+trcData.PSIS_L(cropInd1:cropInd2,:))/2);
    mid = (front+back)/2;
    left = mean(trcData.ASIS_L(cropInd1:cropInd2,:));
    right = mean(trcData.ASIS_R(cropInd1:cropInd2,:));

    %     P = [front(1:2) 0];
    %     Q = [mid(1:2) 0];
    %     R = [];

    %This is rotation from world to body orieantation in mocap
    [~,R] = points2rot([front(1:2) 0],[mid(1:2) 0],[left(1:2) 0]);
    %R = eye(3);
    %     [x y z] = dcm2angle(R,'XYZ');

    %R = eye(3);
    lasis = mean(R*trcData.ASIS_L(cropInd1:cropInd2,:)',2);
    rasis = mean(R*trcData.ASIS_R(cropInd1:cropInd2,:)',2);
    lback = mean(R*trcData.PSIS_L(cropInd1:cropInd2,:)',2);
    rback = mean(R*trcData.PSIS_R(cropInd1:cropInd2,:)',2);
    lankle = (mean(R*trcData.ANKLE_L_LAT(cropInd1:cropInd2,:)',2) + mean(R*trcData.ANKLE_L_MED(cropInd1:cropInd2,:)',2))/2;
    rankle = (mean(R*trcData.ANKLE_R_LAT(cropInd1:cropInd2,:)',2) + mean(R*trcData.ANKLE_R_MED(cropInd1:cropInd2,:)',2))/2;

    %Distance between hips as calculated in

    %@ARTICLE{harrington2007prediction,
    %    author = {Harrington, ME and Zavatsky, AB and Lawson, SEM and Yuan, Z and Theologis,TN},
    %    title = {Prediction of the hip joint centre in adults, children, and patients
    %    with cerebral palsy based on magnetic resonance imaging},
    %    journal = {Journal of biomechanics},
    %    year = {2007},
    %    volume = {40}
    %}

    %Pelvis Width calculated using x and y only, converting back to [mm]
    PW = norm(lasis(1:2) - rasis(1:2))*1000;
    %Pelvis Depth calculated using x and y only
    PD = norm(abs((lasis(1:2)+rasis(1:2))/2 - (lback(1:2)+rback(1:2))/2))*1000;
    %Leg Length
    LL = (norm(rasis-rankle) + norm(lasis-lankle))/2*1000;

    %So going from the middle the joint center is predicted in mm as
    x = -0.24*PD-9.9;
    y = -0.16*PW-0.04*LL-7.1;
    z = 0.28*PD+0.16*PW+7.9;
    % Where x is from middle to front, z is from middle to side, and y is
    % up
    
    %Middle to RIGHT/LEFT hip JOINT CENTERS, converting back to [m]
%     mid2RHC = ([x -z y]/1000)';
%     mid2LHC = ([x z y]/1000)';
%     mid2RHC = ([z x y]/1000)';
%     mid2LHC = ([-z x y]/1000)';

    asismid = (lasis + rasis)/2;
    hipcentre = (lasis + rasis + lback + rback) / 4;
    
    vecRight_Harrington = asismid + rotx(pi/2)*[x; y; z]/1000; % this rotation is to convert harrington's frame into the frame where the calculations have been done
    vecLeft_Harrington = asismid + rotx(pi/2)*[x; y; -z]/1000;
    
        % rotating from Harrington/ISB frame to RL frame
% % %     rot = rotx(-pi/2)*rotz(-pi/2);
% % %     la = [1 2 3];
% % %     rot'*la'
%     mid2RHC = [z; x; y]/1000; % need to rotate in  (rotx -> rotz)
%     mid2LHC = [-z; x; y]/1000;

%     vecR = [x; y; z]/1000;
%     vecL = [x; y; -z]/1000;
%     
%     mid2RHC = R'*vecR;
%     mid2LHC = R'*vecL;
    
    if 0
        figure(4); clf; hold on; grid on;
        plot3(lasis(1),lasis(2),lasis(3),'c*');
        plot3(rasis(1),rasis(2),rasis(3),'c*');
        plot3(lback(1),lback(2),lback(3),'r*');
        plot3(rback(1),rback(2),rback(3),'r*');
        plot3(asismid(1),asismid(2),asismid(3),'k*');
        plot3(hipcentre(1),hipcentre(2),hipcentre(3),'k*');
        plot3(lankle(1),lankle(2),lankle(3),'b*');
        plot3(rankle(1),rankle(2),rankle(3),'b*');
        plot3(vecRight_Harrington(1),vecRight_Harrington(2),vecRight_Harrington(3),'go');
        plot3(vecLeft_Harrington(1),vecLeft_Harrington(2),vecLeft_Harrington(3),'go');
        
        text(lasis(1),lasis(2),lasis(3),'ASIS L');
        text(rasis(1),rasis(2),rasis(3),'ASIS R');
        text(lback(1),lback(2),lback(3),'PSIS L');
        text(rback(1),rback(2),rback(3),'PSIS R');
        text(asismid(1),asismid(2),asismid(3),'');
        text(hipcentre(1),hipcentre(2),hipcentre(3),'HIP BASE');
        text(lankle(1),lankle(2),lankle(3),'ANKLE L');
        text(rankle(1),rankle(2),rankle(3),'ANKLE R');
        text(vecRight_Harrington(1),vecRight_Harrington(2),vecRight_Harrington(3),'HJC R');
        text(vecLeft_Harrington(1),vecLeft_Harrington(2),vecLeft_Harrington(3),'HJC L');
        
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        view([-70 70])
        pbaspect([1 1 1])
    end
    
    % now figure out what the HJC is from the original frame
    lasis = mean(trcData.ASIS_L(cropInd1:cropInd2,:)',2);
    rasis = mean(trcData.ASIS_R(cropInd1:cropInd2,:)',2);
    lback = mean(trcData.PSIS_L(cropInd1:cropInd2,:)',2);
    rback = mean(trcData.PSIS_R(cropInd1:cropInd2,:)',2);
    lankle = (mean(trcData.ANKLE_L_LAT(cropInd1:cropInd2,:)',2) + mean(trcData.ANKLE_L_MED(cropInd1:cropInd2,:)',2))/2;
    rankle = (mean(trcData.ANKLE_R_LAT(cropInd1:cropInd2,:)',2) + mean(trcData.ANKLE_R_MED(cropInd1:cropInd2,:)',2))/2;

    asismid = (lasis + rasis)/2;
    hipcentre = (lasis + rasis + lback + rback) / 4;
    
    vecR = [x; y; z]/1000;
    vecL = [x; y; -z]/1000;
    
    Rint = rotx(-pi/2)*rotz(-pi/2); % transform from Harrington frame to RL frame
    Rin2 = rotz(pi/2); % point2rot didn't actually bring it into Harrington frame, there's still this offset
    mid2RHC = Rin2'*R'*Rint'*vecR;
    mid2LHC = Rin2'*R'*Rint'*vecL;
    
    global2RHC = asismid + mid2RHC; % this rotation is to convert harrington's frame into the frame where the calculations have been done
    global2LHC = asismid + mid2LHC;
    
    if 0
        figure(5); clf; hold on; grid on;
        plot3(lasis(1),lasis(2),lasis(3),'c*');
        plot3(rasis(1),rasis(2),rasis(3),'c*');
        plot3(lback(1),lback(2),lback(3),'r*');
        plot3(rback(1),rback(2),rback(3),'r*');
        plot3(asismid(1),asismid(2),asismid(3),'k*');
        plot3(hipcentre(1),hipcentre(2),hipcentre(3),'k*');
        plot3(lankle(1),lankle(2),lankle(3),'b*');
        plot3(rankle(1),rankle(2),rankle(3),'b*');
        plot3(global2RHC(1),global2RHC(2),global2RHC(3),'go');
        plot3(global2LHC(1),global2LHC(2),global2LHC(3),'go');
        
        text(lasis(1),lasis(2),lasis(3),'ASIS L');
        text(rasis(1),rasis(2),rasis(3),'ASIS R');
        text(lback(1),lback(2),lback(3),'PSIS L');
        text(rback(1),rback(2),rback(3),'PSIS R');
        text(asismid(1),asismid(2),asismid(3),'');
        text(hipcentre(1),hipcentre(2),hipcentre(3),'HIP BASE');
        text(lankle(1),lankle(2),lankle(3),'ANKLE L');
        text(rankle(1),rankle(2),rankle(3),'ANKLE R');
        text(global2RHC(1),global2RHC(2),global2RHC(3),'HJC R');
        text(global2LHC(1),global2LHC(2),global2LHC(3),'HJC L');
        
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        view([33 33])
        pbaspect([1 1 1])
    end

end