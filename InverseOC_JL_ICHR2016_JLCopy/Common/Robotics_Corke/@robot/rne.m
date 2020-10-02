%RNE Compute inverse dynamics via recursive Newton-Euler formulation
%
%	TAU = RNE(ROBOT, Q, QD, QDD)
%	TAU = RNE(ROBOT, [Q QD QDD])
%
% Returns the joint torque required to achieve the specified joint position,
% velocity and acceleration state.
%
% Gravity vector is an attribute of the robot object but this may be 
% overriden by providing a gravity acceleration vector [gx gy gz].
%
%	TAU = RNE(ROBOT, Q, QD, QDD, GRAV)
%	TAU = RNE(ROBOT, [Q QD QDD], GRAV)
%
% An external force/moment acting on the end of the manipulator may also be
% specified by a 6-element vector [Fx Fy Fz Mx My Mz].
%
%	TAU = RNE(ROBOT, Q, QD, QDD, GRAV, FEXT)
%	TAU = RNE(ROBOT, [Q QD QDD], GRAV, FEXT)
%
% where Q, QD and QDD are row vectors of the manipulator state; pos, vel, 
% and accel.
%
% The torque computed also contains a contribution due to armature
% inertia.
%
% RNE can be either an M-file or a MEX-file.  See the manual for details on
% how to configure the MEX-file.  The M-file is a wrapper which calls either
% RNE_DH or RNE_MDH depending on the kinematic conventions used by the robot
% object.
%
% See also: ROBOT, ACCEL, GRAVLOAD, INERTIA.

%
% verified against MAPLE code, which is verified by examples
%

% Copyright (C) 1992-2008, by Peter I. Corke
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


function tau = rne(robot, varargin)
	if robot.mdh == 0,
		tau = rne_dh(robot, varargin{:});
	else
		tau = rne_mdh(robot, varargin{:});
	end
