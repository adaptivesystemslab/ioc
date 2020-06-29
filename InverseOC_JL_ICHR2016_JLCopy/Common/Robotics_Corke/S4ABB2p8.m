% S4ABB2p8  Load kinematic data for an ABB S4 2.8robot 
%
%	S4ABB2p8
%
% Defines the object 'R' in the current workspace which describes the 
% kinematic characterstics of an ABB S4 2.8 robot
% using standard DH conventions.
%
% Also define the vector q0 which corresponds to mastering position.
%
%Wynand Swart
%Mega Robots CC, P/O Box 8412, Pretoria, 0001, South Africa
%Cell: 073-1555-430
%wynand.swart@gmail.com
%30 Sep 2007
%S4 ABB 2.8 robot

L1 = link([ -pi/2   0.188      0        0.9             0]);
L2 = link([ 0       0.95       0        0               0]);
L3 = link([ -pi/2   0.225      0        0               0]);
L4 = link([ pi/2    0          0        1.705           0]);
L5 = link([ -pi/2   0          0        0               0]);
L6 = link([ -pi/2   0          0        0.2             0]);
%##########################################################
%Pose 0; At SYNCHRONISATION position
%##########################################################
q0 = [0     -pi/2         0       0      0     -pi/2];
R=robot({L1 L2 L3 L4 L5 L6});
R.name='S4 ABB 2.8';
%##########################################################
