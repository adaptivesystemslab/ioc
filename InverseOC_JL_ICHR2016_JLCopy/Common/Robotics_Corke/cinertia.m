%CINERTIA Compute the Cartesian (operational space) manipulator inertia matrix
%
%	M = CINERTIA(ROBOT, Q)
%
% Return the n x n inertia matrix which relates Cartesian force/torque to 
% Cartesian acceleration.
% ROBOT is an n-axis robot object and describes the manipulator dynamics and 
% kinematics, and Q is an n element vector of joint state.
%
% See also: INERTIA, ROBOT, RNE.

% MOD HISTORY
% 	4/99 add object support
% $Log: not supported by cvs2svn $
% $Revision: 1.2 $

% Copyright (C) 1993-2008, by Peter I. Corke
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

function Mx = cinertia(robot, q)
	J = jacob0(robot, q);
	Ji = inv(J);
	M = inertia(robot, q);
	Mx = Ji' * M * Ji;
