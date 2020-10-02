% add Corke's toolbox
cd('D:\aslab_gitlab\kalmanfilter\ik_framework\toolboxes\Robotics_Corke')
startup_rvc;

% Model definition with Peter Corke's toolbox
% clear;close all;

% m1 = 1; 
% m2 = 1;
% r1 = [0.005 -.5 0.005];
% r2 = [0.005 -.5 0.005];
% i1 = [1/12*(1*1+0.01*0.01) 0 0;...
%     0 1/12*(0.01*0.01 + 0.01*0.01) 0; 0 0 1/12*(1*1 + 0.05*0.05 )];
% i2 = [1/12*(1*1+0.01*0.01) 0 0;...
%     0 1/12*(0.01*0.01 + 0.01*0.01) 0; 0 0 1/12*(1*1 + 0.01*0.01 )];

l1 = norm([0.3240 0 0]);
m1 = 1.2540;
r1 = [0.1471 -0.0237 -0.0091];
i1 = [    0.0038   -0.0001   -0.0026
   -0.0001    0.0143   -0.0003
   -0.0026   -0.0003    0.0143];

l2 = norm([0.2192 0 0]);
m2 = 0.7410;
r2 = [0.0901 0.0046 0.0042];
i2 = [    0.0007   -0.0004    0.0006
   -0.0004    0.0024    0.0001
    0.0006    0.0001    0.0022];

L(1) = RevoluteMDH('d', 0, 'a', 0, 'alpha', 0, 'offset', 0,...
    'm', 0, 'r', [0 0 0], 'I', zeros(3,3), 'qlim', [-pi*2/3 pi/2]);
L(2) = RevoluteMDH('d', 0, 'a', 0, 'alpha', -pi/2, 'offset', pi/2, ...
    'm', 0, 'r',[0 0 0] , 'I', zeros(3,3), 'qlim', [-pi/2 pi/2]); 
L(3) = RevoluteMDH('d', 11, 'a', 0, 'alpha', pi/2, 'offset', 0,...
    'm', m1, 'r', r1, 'I', i1, 'qlim', [-2*pi/3 pi]); 

L(4) = RevoluteMDH('d', 0, 'a', 0, 'alpha', -pi/2, 'offset', 0,...
    'm', 0, 'r', [0 0 0], 'I', zeros(3, 3), 'qlim', [0 5*pi/6]); 
L(5) = RevoluteMDH('d', 12, 'a', 0, 'alpha', pi/2, 'offset', 0,...
    'm', m2, 'r', r2, 'I', i2, 'qlim', [0 0]);

twolink = SerialLink(L, 'name', 'twolink');
twolink.plot([0 0 0 0 0]);

% Torque with arm up and down to the side
disp(twolink.gravload([0 0 0 0 0]));
disp(twolink.gravload([0 pi/2 0 0 0]));