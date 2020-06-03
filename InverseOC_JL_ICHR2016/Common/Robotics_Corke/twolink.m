%TWOLINK Load kinematic and dynamic data for a simple 2-link mechanism
%
%	TWOLINK
%
% Defines the object 'tl' in the current workspace which describes the 
% kinematic and dynamic characterstics of a simple planar 2-link mechanism.
%
% Example based on Fig 3-6 (p73) of Spong and Vidyasagar (1st edition).  
% It is a planar mechanism operating in the XY (horizontal) plane and is 
% therefore not affected by gravity.
%
% Assume unit length links with all mass (unity) concentrated at the joints.
%
% Also define the vector qz = [0 0] which corresponds to the zero joint
% angle configuration.
%
% See also: DH, DYN, PUMA560, PUMA560AKB, STANFORD.

% Copyright (C) 2000-2008, by Peter I. Corke
%
% This file is part of The Robotics Toolbox for Matlab (RTB).
% 
% RTB is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% RTB is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with RTB.  If not, see <http://www.gnu.org/licenses/>.

twolink_dh = [
% alpha A	theta	D	sigma	m	rx	ry	rz	Ixx	Iyy	Izz	Ixy	Iyz	Ixz	Jm	G
  0     1         0     0         0     1       1       0       0       0       0       0       0       0       0        0      1
  0     1         0     0         0     1       1       0       0       0       0       0       0       0       0        0      1
];

tl = robot(twolink_dh);
tl.name = 'Simple two link';
qz = [0 0];
