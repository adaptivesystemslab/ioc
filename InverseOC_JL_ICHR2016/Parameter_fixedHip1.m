%% Geometrical parameters
param.NF=9;param.NJ=9;


% Mass Parameters
param.M1=  (30.85/100)*param.body_weight; % L5: [Dumas2007: LJC to projection of HJC in sagittal plane + 1/2 torso] hip to shoulder
param.M2=2*(12.3/100)*param.body_weight; % L3: [Dumas2007: HJC to KJC, hip to knee] 
param.M3=2*(4.8/100)*param.body_weight;  % L4: [Dumas2007: KJC to AJC, knee to ankle ]
param.M4=2*(1.2/100)*param.body_weight*5;  % L5: [Dumas2007: CAL to TTII, ie ankle to toe]
param.M5=0;
param.M6=0;
param.M7=0;
param.M8=0;
param.M9=0;
param.M10=0;

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