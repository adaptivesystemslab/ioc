%INERTIA Compute the manipulator inertia matrix
%
%	INERTIA(ROBOT, Q)
%
% Returns the n x n symmetric inertia matrix which relates joint torque 
% to joint acceleration.
% ROBOT describes the manipulator dynamics and kinematics, and Q is
% an n element vector of joint state.
%
% See also: RNE, CINERTIA, ITORQUE, CORIOLIS, GRAVLOAD.

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

function M = inertia(robot, q)
	n = robot.n;

	if numel(q) == robot.n,
		q = q(:)';
	end

	M = zeros(n,n,0);
	for Q = q',
		m = rne(robot, ones(n,1)*Q', zeros(n,n), eye(n), [0;0;0]);
		M = cat(3, M, m);
	end
