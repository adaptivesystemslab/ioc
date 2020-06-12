clearvars;

m1 = 1.2540;
l1 = [0.3240  0       0];
r1 = [0.1471 -0.0237 -0.0091];
i1 = [    0.0038   -0.0001   -0.0026
   -0.0001    0.0143   -0.0003
   -0.0026   -0.0003    0.0143];

l1f = l1(1);
r1f = [r1(1) 0 0];

%% Load Corke Model
% L(1) = RevoluteMDH('d', 0, 'a', 0, 'alpha', 0, 'offset', 0,...
%     'm', 0, 'r', [0 0 0], 'I', zeros(3,3), 'qlim', [-pi*2/3 pi/2]);
% L(2) = RevoluteMDH('d', 0, 'a', 0, 'alpha', -pi/2, 'offset', pi/2, ...
%     'm', 0, 'r',[0 0 0] , 'I', zeros(3,3), 'qlim', [-pi/2 pi/2]); 
% L(3) = RevoluteMDH('d', l1f, 'a', 0, 'alpha', pi/2, 'offset', 0,...
%     'm', m1, 'r', r1f, 'I', i1, 'qlim', [-2*pi/3 pi]); 
% L(4) = RevoluteMDH('d', 0, 'a', 0, 'alpha', -pi/2, 'offset', 0,...
%     'm', 0, 'r', [0 0 0], 'I', zeros(3, 3), 'qlim', [0 5*pi/6]); 

L(1) = RevoluteMDH('d', 0, 'a', 0, 'alpha', pi/2, 'offset', pi/2,...
    'm', m1, 'r', r1f, 'I', i1, 'qlim', [-2*pi/3 pi]); 

mdl = SerialLink(L, 'name', 'mdl');

q = [0];
dq = [0];
ddq = [0];

mdl.plot(q);

tau_corke = mdl.rne(q, dq, ddq)
tau_hand = m1*r1f*9.81