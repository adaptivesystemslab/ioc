% param.body_weight=75;
% param.body_height=1.75;
%% Geometrical parameters
% param.L1=param.body_height*(1/1.77)*0.219; % X offset from ankle
% param.L2=param.body_height*(1/1.77)*0.06; % Y offset from ankle
% param.L3=param.body_height*(1/1.77)*0.433; % ankle to knee
% param.L4=param.body_height*(1/1.77)*0.432; % knee to hip
% param.L5=param.body_height*(1/1.77)*0.094; % hip to pelvis (from ASIS/PSIS)
% param.L6=param.body_height*(1/1.77)*0.151;  % pelvis to lower torso (from T10/STR)
% param.L7=param.body_height*(1/1.77)*0.334;  % mid torso to higher torso (C7/CLAV)
% param.L8=param.body_height*(1/1.77)*0.271; % higher torso to elbow    
param.L9=param.body_height*(1/1.77)*0.283;%forearm
param.L10=param.body_height*(1/1.77)*0.189;%hand
param.L11=param.body_height*(1/1.77)*0.277;%Head&Neck

param.NF=9;param.NJ=9;
param.CALtoTT=0.265; % originally CALtoTT (L0)

% param.L1=abs(Markers.RANK(1,1)); % RANK, frame 1, X
% param.L2=abs(Markers.RANK(1,2)); % RANK, frame 1, Y
% param.L3=norm(Markers.RANK(1,1:2)-Markers.RKNE(1,1:2)); % RANK-RKNE, frame 1, X&Y
% param.L4=norm(Markers.RKNE(1,1:2)-Markers.RHJC(1,1:2));
% param.L5=norm(Markers.RHJC(1,1:2)-Markers.Mid_Pelvis(1,1:2)); % hip to pelvis (from ASIS/PSIS)

% Mass Parameters
param.M1=2*(1.2/100)*param.body_weight;  % L1: [Dumas2007: CAL to TTII, ie ankle to toe]
param.M2=2*(4.8/100)*param.body_weight;  % L3: [Dumas2007: KJC to AJC, knee to ankle ]
param.M3=2*(12.3/100)*param.body_weight; % L4: [Dumas2007: HJC to KJC, hip to knee] 
param.M4=(14.2/100)*param.body_weight;   % L5: [Dumas2007: LJC to projection of HJC in sagittal plane] hip to pelvis
param.M5=(2.9/100)*param.body_weight; % 
param.M6=(30.4/100)*param.body_weight;
param.M7=2*(2.4/100)*param.body_weight;
param.M8=2*(1.7/100)*param.body_weight;
param.M9=2*(0.6/100)*param.body_weight;
param.M10=(6.7/100)*param.body_weight;

param.Mtot=param.M1+param.M2+param.M3+param.M4+param.M5+param.M6+param.M7+param.M8+param.M9+param.M10;



% Center of Mass
R1=rotz(90);

param.COM1=param.CALtoTT*[-43.6;-2.5;0]/100;
param.COM2=R1*param.L3*[-4.8;-41.0;0]/100;
param.COM3=R1*param.L4*[-4.1;-42.9;0]/100;
param.COM4=R1*param.L5*[2.8;-28.0;0]/100;
param.COM5=R1*param.L6*[17.3;-36.1;0]/100;
param.COM6=R1*param.L7*[0;-55.5;0]/100;
param.COM7=-R1*param.L8*[1.7;-45.2;0]/100;
% param.COM8=(param.M1*param.COM1+param.M2*param.COM2+param.M3*param.COM3+param.M4*param.COM4+param.M5*param.COM5+param.M6*param.COM6+param.M7*param.COM7)/param.Mtot;

 for i=1:7
    eval([ 'param.MX' num2str(i) '=param.M' num2str(i) '*param.COM' num2str(i) '(1,1);']);
    eval([ 'param.MY' num2str(i) '=param.M' num2str(i) '*param.COM' num2str(i) '(2,1);']);
 end
    



%% Inertia Parameters


    %From COM to origin in Dumas frame 

IX_coord=[param.CALtoTT*43.6;...
          param.L3*-4.8;...
          param.L4*-4.1;...
          param.L5*-33.6;...
          param.L6*17.3;...
          param.L7*0;...
          param.L8*1.7;...
          param.L9*1.0;...
          param.L10*3.5;...
          param.L11*2]/100;
IY_coord=[param.CALtoTT*-2.5;...
          param.L3*-41.0;...
          param.L4*-42.9;...
          param.L5*-14.9;...
          param.L6*-36.1;...
          param.L7*-55.5;...
          param.L8*-45.2;...
          param.L9*-41.7;...
          param.L10*-35.7;...
          param.L11*53.6]/100;
IZ_coord=[param.CALtoTT*-0.7;...
          param.L3*0.7;...
          param.L4*3.3;...
          param.L5*3.2;...
          param.L6*-3.3;...
          param.L7*-0.4;...
          param.L8*-2.6;...
          param.L9*1.4;...
          param.L10*3.2;...
          param.L11*0.1]/100;
%%
for i=1:10
    eval(['IMat_huygens' num2str(i) '=param.M' num2str(i) '*[IY_coord(' num2str(i) ',1)^2+ IZ_coord(' num2str(i) ',1)^2 , -IX_coord(' num2str(i) ',1)*IY_coord(' num2str(i) ',1),-IX_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1);-IX_coord(' num2str(i) ',1)*IY_coord(' num2str(i) ',1),IX_coord(' num2str(i) ',1)^2+IZ_coord(' num2str(i) ',1)^2,-IY_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1);-IX_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1),-IY_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1),IX_coord(' num2str(i) ',1)^2+IY_coord(' num2str(i) ',1)^2];'])
end
%%
link_1=(IMat_huygens1+param.M1*[(0.11*(param.CALtoTT))^2 (0.09*(param.CALtoTT))^2 (0.06*1i*(param.CALtoTT))^2;...
                    (0.09*(param.CALtoTT))^2 (0.25*(param.CALtoTT))^2 0;...
                    (0.06*1i*(param.CALtoTT))^2 0 (0.25*(param.CALtoTT))^2]);
                
link_2=(IMat_huygens2+param.M2*[(0.28*param.L3)^2 (0.04*1i*param.L3)^2 (0.02*1i*param.L3)^2;...
                    (0.04*1i*param.L3)^2 (0.10*param.L3)^2 (0.05*param.L3)^2;...
                    (0.02*1i*param.L3)^2 (0.05*param.L3)^2 (0.28*param.L3)^2]);

                
link_3=(IMat_huygens3+param.M3*[(0.29*param.L4)^2 (0.07*param.L4)^2 (0.02*1i*param.L4)^2;...
                    (0.07*param.L4)^2 (0.15*param.L4)^2 (0.07*1i*param.L4)^2;...
                    (0.02*1i*param.L4)^2 (0.07*1i*param.L4)^2 (0.30*param.L4)^2]);

link_4=(IMat_huygens4+param.M4*[(0.42*param.L5)^2 (0.10*1i*param.L5)^2 (0.05*1i*param.L5)^2;...
                    (0.10*1i*param.L5)^2 (0.44*param.L5)^2 (0.03*1i*param.L5)^2;...
                    (0.05*1i*param.L5)^2 (0.03*1i*param.L5)^2 (0.40*param.L5)^2]);
                
link_5=(IMat_huygens5+param.M5*[(0.54*param.L6)^2 (0.11*param.L6)^2 (0.06*1i*param.L6)^2;...
                    (0.11*param.L6)^2 (0.66*param.L6)^2 (0.05*1i*param.L6)^2;...
                    (0.06*1i*param.L6)^2 (0.05*1i*param.L6)^2 (0.40*param.L6)^2]);

link_6=(IMat_huygens6+param.M6*[(0.42*param.L7)^2 (0.11*1i*param.L7)^2 (0.01*param.L7)^2;...
                    (0.11*1i*param.L7)^2 (0.33*param.L7)^2 (0.03*param.L7)^2;...
                    (0.01*param.L7)^2 (0.03*param.L7)^2 (0.36*param.L7)^2]);                
                
                
link_7=(IMat_huygens7+param.M7*[(0.31*param.L8)^2 (0.06*param.L8)^2 (0.05*param.L8)^2;...
                    (0.06*param.L8)^2 (0.14*param.L8)^2 (0.02*param.L8)^2;...
                    (0.05*param.L8)^2 (0.02*param.L8)^2 (0.32*param.L8)^2]); 
                

link_8=(IMat_huygens8+param.M8*[(0.28*param.L9)^2 (0.03*param.L9)^2 (0.02*param.L9)^2;...
                                (0.03*param.L9)^2 (0.11*param.L9)^2 (0.08*1i*param.L9)^2;...
                                (0.02*param.L9)^2 (0.08*1i*param.L9)^2 (0.27*param.L9)^2]);    
                            
                            
link_9=(IMat_huygens9+param.M9*[(0.26*param.L10)^2 (0.09*param.L10)^2 (0.07*param.L10)^2;...
                                (0.09*param.L10)^2 (0.16*param.L10)^2 (0.08*1i*param.L10)^2;...
                                (0.07*param.L10)^2 (0.08*1i*param.L10)^2 (0.24*param.L10)^2]);

link_10=(IMat_huygens10+param.M10*[(0.28*param.L11)^2 (0.07*1i*param.L11)^2 (0.02*1i*param.L11)^2;...
                                (0.07*1i*param.L11)^2 (0.21*param.L11)^2 (0.03*param.L11)^2;...
                                (0.02*1i*param.L11)^2 (0.03*param.L11)^2 (0.30*param.L11)^2]);

% From origin to the right point            
                
    IX_coord=[-param.L1;0;0;0;0;0;0;0;0;0];
    IY_coord=[param.L2;param.L3;param.L4;param.L5;param.L6;param.L7;0;-param.L8;-(param.L8+param.L9);param.L7];
    IZ_coord=[0;0;0;0;0;0;0;0;0;0];

for i=1:10
    eval(['IMat_huygens' num2str(i) '=param.M' num2str(i) '*[IY_coord(' num2str(i) ',1)^2+ IZ_coord(' num2str(i) ',1)^2 , -IX_coord(' num2str(i) ',1)*IY_coord(' num2str(i) ',1),-IX_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1);-IX_coord(' num2str(i) ',1)*IY_coord(' num2str(i) ',1),IX_coord(' num2str(i) ',1)^2+IZ_coord(' num2str(i) ',1)^2,-IY_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1);-IX_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1),-IY_coord(' num2str(i) ',1)*IZ_coord(' num2str(i) ',1),IX_coord(' num2str(i) ',1)^2+IY_coord(' num2str(i) ',1)^2];'])
    eval(['link_' num2str(i) '=link_' num2str(i) '+IMat_huygens' num2str(i) ';']);
end

% Forearms and hands on the frame of the shoulder
link_7=link_7+link_8+link_9;
% Head on the frame of the upper torso
link_6=link_6+link_10;

%projection matrix in Model frame    
    %rotation matrix between Dumas'frame and Model's frame 
    R1=rotz(90);
for i=1:7
    eval(['link_' num2str(i) '=transpose(R1)*link_' num2str(i) '*R1;'])
end

                
 for i=1:7
    eval([ 'param.ZZ' num2str(i) '=link_' num2str(i) '(3,3);']);
 end
    
% Update of the total mass for the arm+forearm+hand
param.M7=param.M7+param.M8+param.M9;
% Update of the total mass for the upper trunk+Head&Neck
param.M6=param.M6+param.M10;

% save('param_structure','param')


gravityVal = -9.81;