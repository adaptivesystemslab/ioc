% param.body_weight=75;
% param.body_height=1.75;
%% Geometrical parameters
% param.L1=param.body_height*(1/1.77)*0.277;%Head&Neck
% param.L2=param.body_height*(1/1.77)*0.189;%hand
% param.L3=param.body_height*(1/1.77)*0.283;%forearm
% param.L4=param.body_height*(1/1.77)*0.271; % higher torso to elbow
% param.L5=param.body_height*(1/1.77)*0.334;  % mid torso to higher torso (C7/CLAV)
% param.L6=param.body_height*(1/1.77)*0.151; % pelvis to lower torso (from T10/STR)
% param.L7=param.body_height*(1/1.77)*0.094; % hip to pelvis (from ASIS/PSIS)

param.L1 = 0; % shoulder X
param.L2 = 0; % shoulder Y


param.NF=9;param.NJ=9;
 % originally CALtoTT (L12)
  
% q2: ankle to knee joint
% q3: knee to hip joint
% q4: hip to torso joint

% param.L1=abs(Markers.RANK(1,1)); % RANK, frame 1, X
% param.L2=abs(Markers.RANK(1,2)); % RANK, frame 1, Y
% param.L3=norm(Markers.RANK(1,1:2)-Markers.RKNE(1,1:2)); % RANK-RKNE, frame 1, X&Y
% param.L4=norm(Markers.RKNE(1,1:2)-Markers.RHJC(1,1:2));
% param.L5=norm(Markers.RHJC(1,1:2)-Markers.Mid_Pelvis(1,1:2)); % hip to pelvis (from ASIS/PSIS)

% Mass Parameters
% param.M1=2*(1.2/100)*param.body_weight;  % L1: [Dumas2007: CAL to TTII, ie ankle to toe]
% param.M2=2*(4.8/100)*param.body_weight;  % L3: [Dumas2007: KJC to AJC, knee to ankle ]
% param.M3=2*(12.3/100)*param.body_weight; % L4: [Dumas2007: HJC to KJC, hip to knee] 
% param.M4=(14.2/100)*param.body_weight;   % L5: [Dumas2007: LJC to projection of HJC in sagittal plane] hip to pelvis

% marker{1} = (marker_shoulder - markerOffset)/1000;
% marker{2} = (marker_hip - markerOffset)/1000;
% marker{3} = (marker_knee - markerOffset)/1000;
% marker{4} = (marker_ankle - markerOffset)/1000;
% marker{5} = (marker_toe - markerOffset)/1000;
%   param.L1=abs(Markers{2}(1,1)); % shoulder x
%   param.L2=abs(Markers{2}(1,2)); % shoulder y
%   param.L3=norm(Markers{2}(1,1:2)-Markers{3}(1,1:2)); % hip to knee
%   param.L4=norm(Markers{3}(1,1:2)-Markers{4}(1,1:2)); % knee to ankle 
%   param.L5=norm(Markers{4}(1,1:2)-Markers{5}(1,1:2)); % ankle to toe

% Mass Parameters
param.M1=  (14.2/100)*param.body_weight; % L1: [Dumas2007: LJC to projection of HJC in sagittal plane] hip to pelvis
param.M2=2*(12.3/100)*param.body_weight; % L3: [Dumas2007: HJC to KJC, hip to knee] 
param.M3=2*(4.8/100)*param.body_weight;  % L4: [Dumas2007: KJC to AJC, knee to ankle ]
param.M4=2*(1.2/100)*param.body_weight;  % L5: [Dumas2007: CAL to TTII, ie ankle to toe]
param.M5=2*(1.2/100)*param.body_weight;
param.M6=2*(1.2/100)*param.body_weight;
param.M7=2*(1.2/100)*param.body_weight;
param.M8=2*(1.2/100)*param.body_weight;
param.M9=2*(1.2/100)*param.body_weight;
param.M10=2*(1.2/100)*param.body_weight;

param.Mtot=param.M1+param.M2+param.M3+param.M4+param.M5+param.M6+param.M7+param.M8+param.M9+param.M10;
 
% Center of Mass
R1=rotz(-90);
  
param.COM1=-R1*param.L1*[0;0;0]/100; 
param.COM2= R1*param.L3*-[-4.1 ;-42.9;0]/100;
param.COM3= R1*param.L4*-[-4.8 ;-41.0;0]/100; 
param.COM4= R1*param.L5*-[ 43.6; -2.5;0]/100;
param.COM5= R1*param.L7*[0;0;0]/100;
param.COM6= R1*param.L8*[0;0;0]/100;
param.COM7=    param.L9*[0;0;0]/100; 

% param.COM8=(param.M1*param.COM1+param.M2*param.COM2+param.M3*param.COM3+param.M4*param.COM4+param.M5*param.COM5+param.M6*param.COM6+param.M7*param.COM7)/param.Mtot;

 for i=1:7
    eval([ 'param.MX' num2str(i) '=param.M' num2str(i) '*param.COM' num2str(i) '(1,1);']);
    eval([ 'param.MY' num2str(i) '=param.M' num2str(i) '*param.COM' num2str(i) '(2,1);']);
 end
    



%% Inertia Parameters


    %From COM to origin in Dumas frame 

IX_coord=[  0;...
            0;...
            param.L3*-4.1;...
            param.L4*-4.8;...
            param.L5*43.6;...
            0;...
            0;...
            0;...
            0;...
            0;]/100;
IY_coord=[  0;...
            0;...
            param.L3*-42.9;...
            param.L4*-41.0;...
            param.L5*-2.5;...
            0;...
            0;...   
            0;...
            0;...
            0]/100;
IZ_coord=[  0;...
            0;...
            param.L3*3.3;...
            param.L4*0.7;...
            param.L5*-0.7;...
            0;...
            0;...
            0;...
            0;...
            0]/100;
%%
for i=1:10
    eval(['IMat_huygens' num2str(i) '=param.M' num2str(i) '*[IY_coord(' num2str(i) ',1)^2+ IZ_coord(' num2str(i) ',1)^2 , -IX_coord(' num2str(i) ',1)*IY_coord(' num2str(i) ',1),-IX_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1);-IX_coord(' num2str(i) ',1)*IY_coord(' num2str(i) ',1),IX_coord(' num2str(i) ',1)^2+IZ_coord(' num2str(i) ',1)^2,-IY_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1);-IX_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1),-IY_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1),IX_coord(' num2str(i) ',1)^2+IY_coord(' num2str(i) ',1)^2];'])
end
%%
                      
link_1=(IMat_huygens1+param.M1*zeros(3,3));

link_2=(IMat_huygens2+param.M2*[(0.29*param.L3)^2    (0.07*param.L3)^2    (0.02*1i*param.L3)^2;...
                                (0.07*param.L3)^2    (0.15*param.L3)^2    (0.07*1i*param.L3)^2;...
                                (0.02*1i*param.L3)^2 (0.07*1i*param.L3)^2 (0.30*param.L3)^2]);
                            
link_3=(IMat_huygens3+param.M3*[(0.28*param.L4)^2    (0.04*1i*param.L4)^2 (0.02*1i*param.L4)^2;...
                                (0.04*1i*param.L4)^2 (0.10*param.L4)^2    (0.05*param.L4)^2;...
                                (0.02*1i*param.L4)^2 (0.05*param.L4)^2    (0.28*param.L4)^2]); 
                
link_4=(IMat_huygens4+param.M4*[(0.11*(param.L5))^2 (0.09*(param.L5))^2 (0.06*1i*(param.L5))^2;...
                                (0.09*(param.L5))^2 (0.25*(param.L5))^2 0;...
                                (0.06*1i*(param.L5))^2 0                (0.25*(param.L5))^2]);
                            
link_5=(IMat_huygens5+param.M5*zeros(3,3));
link_6=(IMat_huygens6+param.M6*zeros(3,3));
link_7=(IMat_huygens7+param.M7*zeros(3,3));
link_8=(IMat_huygens8+param.M8*zeros(3,3));
link_9=(IMat_huygens9+param.M9*zeros(3,3));
link_10=(IMat_huygens10+param.M10*zeros(3,3));
                                   
% From origin to the right point                 
IX_coord=[param.L1;0;0;0;0;0;0;0;0;0];
IY_coord=[param.L2;param.L3;param.L4;param.L5;0;0;0;0;0;0];
IZ_coord=[0;0;0;0;0;0;0;0;0;0];

for i=1:10
    eval(['IMat_huygens' num2str(i) '=param.M' num2str(i) '*[IY_coord(' num2str(i) ',1)^2+ IZ_coord(' num2str(i) ',1)^2 , -IX_coord(' num2str(i) ',1)*IY_coord(' num2str(i) ',1),-IX_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1);-IX_coord(' num2str(i) ',1)*IY_coord(' num2str(i) ',1),IX_coord(' num2str(i) ',1)^2+IZ_coord(' num2str(i) ',1)^2,-IY_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1);-IX_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1),-IY_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1),IX_coord(' num2str(i) ',1)^2+IY_coord(' num2str(i) ',1)^2];'])
    eval(['link_' num2str(i) '=link_' num2str(i) '+IMat_huygens' num2str(i) ';']);
end


%projection matrix in Model frame    
%rotation matrix between Dumas'frame and Model's frame 
for i=1:7
    eval(['link_' num2str(i) '=transpose(R1)*link_' num2str(i) '*R1;'])
end
                
 for i=1:7
    eval([ 'param.ZZ' num2str(i) '=link_' num2str(i) '(3,3);']);
 end
    
gravityVal = 9.81;