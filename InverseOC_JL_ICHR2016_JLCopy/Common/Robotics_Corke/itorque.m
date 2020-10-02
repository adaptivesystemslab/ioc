%ITORQUE Compute the manipulator inertia torque
%
%	TAUI = ITORQUE(ROBOT, Q, QDD)
%
% Returns the n-element inertia torque vector at the specified pose and 
% acceleration, that is,
% 	TAUI = INERTIA(Q)*QDD
%
% ROBOT describes the manipulator dynamics and kinematics.
% If Q and QDD are row vectors, the result is a row vector of joint torques.
% If Q and QDD are matrices, each row is interpretted as a joint state 
% vector, and the result is a matrix each row being the corresponding joint 
% torques.
% 
% If ROBOT contains non-zero motor inertia then this will included in the
% result.
%
% See also: RNE, CORIOLIS, INERTIA, GRAVLOAD.

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

function it = itorque(robot, q, qdd)
	it = rne(robot, q, zeros(size(q)), qdd, [0;0;0]);
