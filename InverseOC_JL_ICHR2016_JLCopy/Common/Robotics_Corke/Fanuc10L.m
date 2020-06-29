% FANUC10L  Load kinematic data for a Fanuc AM120iB/10L robot 
%
%	Fanuc10L
%
% Defines the object 'R' in the current workspace which describes the 
% kinematic characterstics of a Fanuc AM120iB/10L
% using standard DH conventions.
%
% Also define the vector q0 which corresponds to mastering position.
%
%Wynand Swart
%Mega Robots CC, P/O Box 8412, Pretoria, 0001, South Africa
%Cell: 073-1555-430
%wynand.swart@gmail.com
%30 Sep 2007
%Fanuc AM120iB/10L robot

L1 = link([-pi/2   0.15   0         0       0  ]);
L2 = link([ pi     0.77   0         0       0  ]);
L3 = link([-pi/2   0.1    0         0       0  ]);
L4 = link([ pi/2   0      0        -0.96    0  ]);
L5 = link([-pi/2   0      0         0       0  ]);
L6 = link([ 0      0      0        -0.1     0  ]);
%##########################################################
%Pose 0; At MASTERING position;
%##########################################################
q0 =[0   -pi/2   0   0   0   0];
R=robot({L1 L2 L3 L4 L5 L6});
R.name='Fanuc AM120iB/10L';
%##########################################################
