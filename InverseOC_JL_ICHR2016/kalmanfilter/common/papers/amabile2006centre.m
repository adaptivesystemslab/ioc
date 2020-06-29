function [RSC, LSC] = amabile2006centre(trcData, cropInd1, cropInd2, subjectHeight)
% https://www.sciencedirect.com/science/article/pii/S0021929016303633#bib31
    %Rotation calculate based on the mean of mid->front and right for y
    %Here we get the world to base transform to match Mocap Data
    front = mean((trcData.ASIS_R(cropInd1:cropInd2,:)+trcData.ASIS_L(cropInd1:cropInd2,:))/2);
    back = mean((trcData.PSIS_R(cropInd1:cropInd2,:)+trcData.PSIS_L(cropInd1:cropInd2,:))/2);
    mid = (front+back)/2;
    left = mean(trcData.ASIS_L(cropInd1:cropInd2,:));
    right = mean(trcData.ASIS_R(cropInd1:cropInd2,:));

    %This is rotation from world to body orieantation in mocap
    [~,R] = points2rot([front(1:2) 0],[mid(1:2) 0],[left(1:2) 0]);
    
%     R = rotx(pi/2)'*R;
    
    C7 = mean(R*trcData.C7(cropInd1:cropInd2,:)',2);
    
    IJ = mean(R*trcData.STERN(cropInd1:cropInd2,:)',2) - C7;  
    PX = mean(R*trcData.DIAPHRAM(cropInd1:cropInd2,:)',2) - C7;
    T8 = mean(R*trcData.SCAPULA(cropInd1:cropInd2,:)',2) - C7;
    ACL = mean(R*trcData.SHOULDER_L(cropInd1:cropInd2,:)',2) - C7;
    ACR = mean(R*trcData.SHOULDER_R(cropInd1:cropInd2,:)',2) - C7;
    ASISR = mean(R*trcData.ASIS_R(cropInd1:cropInd2,:)',2) - C7;
    ASISL = mean(R*trcData.ASIS_L(cropInd1:cropInd2,:)',2) - C7;
    
    C7 = zeros(size(C7));
    
    M3 = mean([IJ PX], 2);
    M4 = mean([C7 T8], 2); 
    D_x = norm(M4 - M3);
    
    MAC = mean([ACL ACR], 2);  
    MH = mean([ASISR ASISL], 2);  
    D_y = norm(MAC - MH);
     
    D_z = 0.129*subjectHeight;
    
    A_coeff = mean([-0.61 -0.61 -0.61 -0.61 -0.63 -0.66]);
    B_coeff = mean([ 0    -0.01  0    -0.01 -0.02 -0.02]);
    C_coeff = mean([ 0.56  0.57  0.56  0.57  0.78  0.78]);
    
    x = -A_coeff*D_x;
    y = B_coeff*D_y;
    z = C_coeff*D_z;

    RSC = C7+[x z y]';
%     LSC = C7+[x y z]';

        vis.addMarker('RSC', RSC, colour3); % vis.addMarker('LSC', LSC, colour3);
        vis.update();
    
    if 0
        vis = rlVisualizer('vis',640,480);
        %             obj.model.position = 0;
        colour1 = [0 0 1 1];
        colour2 = [0 1 0 1];
        colour3 = [1 0 0 1];
        colour4 = [0 0 0 0.2];
        
        vis.addMarker('x-axis', [1 0 0], colour4);
        vis.addMarker('y-axis', [0 1 0], colour4);
        vis.addMarker('z-axis', [0 0 1], colour4);
        
        vis.addMarker('ASIS_R', ASISR, colour1); vis.addMarker('ASIS_L', ASISL, colour1);
        vis.addMarker('IJ', IJ, colour1);        vis.addMarker('PX', PX, colour1);
        vis.addMarker('C7', C7, colour1);        vis.addMarker('T8', T8, colour1);
        vis.addMarker('ACL', ACL, colour1);      vis.addMarker('ACR', ACR, colour1);
        
        vis.addMarker('M3', M3, colour2);        vis.addMarker('M4', M4, colour2);
        vis.addMarker('MAC', MAC, colour2);        vis.addMarker('MH', MH, colour2);
        
%         vis.addMarker('RSC', RSC, colour3); vis.addMarker('LSC', LSC, colour3);

        vis.update();
    end
    
    if 0
        figure; hold on; grid on
        plot3(ASISR(1),ASISR(2),ASISR(3),'c*');
        plot3(ASISL(1),ASISL(2),ASISL(3),'c*');
        text(ASISR(1),ASISR(2),ASISR(3),'ASIS_R');
        text(ASISL(1),ASISL(2),ASISL(3),'ASIS_L');
        
        plot3(ACL(1),ACL(2),ACL(3),'k*');
        plot3(ACR(1),ACR(2),ACR(3),'k*');
        text(ACL(1),ACL(2),ACL(3),'ACL');
        text(ACR(1),ACR(2),ACR(3),'ACR');
        
        plot3(IJ(1),IJ(2),IJ(3),'r*');
        plot3(PX(1),PX(2),PX(3),'r*');
        plot3(C7(1),C7(2),C7(3),'r*');
        plot3(T8(1),T8(2),T8(3),'r*');
        text(IJ(1),IJ(2),IJ(3),'IJ');
        text(PX(1),PX(2),PX(3),'PX');
        text(C7(1),C7(2),C7(3),'C7');
        text(T8(1),T8(2),T8(3),'T8');
        
        plot3(M3(1),M3(2),M3(3),'b*');
        plot3(M4(1),M4(2),M4(3),'b*');
        plot3(MAC(1),MAC(2),MAC(3),'b*');
        plot3(MH(1),MH(2),MH(3),'b*');
        text(M3(1),M3(2),M3(3),'M3');
        text(M4(1),M4(2),M4(3),'M4');
        text(MAC(1),MAC(2),MAC(3),'MAC');
        text(MH(1),MH(2),MH(3),'MH');
        
        x = (1+A_coeff)*D_x;
        y = B_coeff*D_y;
        z = C_coeff*D_z;
        C7
        RSC = C7+[x y z]'
%         ACR
        plot3(RSC(1),RSC(2),RSC(3),'ko');
        
        xlabel('x');
        ylabel('y');
        zlabel('z');
        daspect([1 1 1]);
    end
    


% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     %R = eye(3);
%     lasis = mean(R*trcData.ASIS_L(cropInd1:cropInd2,:)',2);
%     rasis = mean(R*trcData.ASIS_R(cropInd1:cropInd2,:)',2);
%     lback = mean(R*trcData.PSIS_L(cropInd1:cropInd2,:)',2);
%     rback = mean(R*trcData.PSIS_R(cropInd1:cropInd2,:)',2);
%     lankle = (mean(R*trcData.ANKLE_L_LAT(cropInd1:cropInd2,:)',2) + mean(R*trcData.ANKLE_L_MED(cropInd1:cropInd2,:)',2))/2;
%     rankle = (mean(R*trcData.ANKLE_R_LAT(cropInd1:cropInd2,:)',2) + mean(R*trcData.ANKLE_R_MED(cropInd1:cropInd2,:)',2))/2;
% 
%     %     figure(2); clf; hold on; grid on;
%     %     plot3(lasis(1),lasis(2),lasis(3),'c*');
%     %     plot3(rasis(1),rasis(2),rasis(3),'r*');
%     %     plot3(lback(1),lback(2),lback(3),'g*');
%     %     plot3(rback(1),rback(2),rback(3),'b*');
%     %     xlabel('X');
%     %     ylabel('Y');
% 
%     %Distance between hips as calculated in
% 
%     %@ARTICLE{harrington2007prediction,
%     %    author = {Harrington, ME and Zavatsky, AB and Lawson, SEM and Yuan, Z and Theologis,TN},
%     %    title = {Prediction of the hip joint centre in adults, children, and patients
%     %    with cerebral palsy based on magnetic resonance imaging},
%     %    journal = {Journal of biomechanics},
%     %    year = {2007},
%     %    volume = {40}
%     %}
% 
%     %Pelvis Width calculated using x and y only, converting back to [mm]
%     PW = norm(lasis(1:2) - rasis(1:2))*1000;
%     %Pelvis Depth calculated using x and y only
%     PD = norm(abs((lasis(1:2)+rasis(1:2))/2 - (lback(1:2)+rback(1:2))/2))*1000;
%     %Leg Length
%     LL = (norm(rasis-rankle) + norm(lasis-lankle))/2*1000;
% 
%     %So going from the middle the joint center is predicted in mm as
%     x = -0.24*PD-9.9;
%     y = -0.16*PW-0.04*LL-7.1;
%     z = 0.28*PD+0.16*PW+7.9;
%     % Where x is from middle to front, z is from middle to side, and y is
%     % up
% 
%     %Middle to RIGHT/LEFT hip JOINT CENTERS, converting back to [m]
%     mid2RHC = ([z x y]/1000)';
%     mid2LHC = ([-z x y]/1000)';
end