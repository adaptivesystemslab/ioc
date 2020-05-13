%FKINE  Forward robot kinematics for serial link manipulator
%
%	TR = FKINE(ROBOT, Q)
%
% Computes the forward kinematics for each joint space point defined by Q.
% ROBOT is a robot object.
%
% For an n-axis manipulator Q is an n element vector or an m x n matrix of
% robot joint coordinates.
% 
% If Q is a vector it is interpretted as the generalized joint coordinates, and
% FKINE(ROBOT, Q) returns a 4x4 homogeneous transformation for the tool of
% the manipulator.
%
% If Q is a matrix, the rows are interpretted as the generalized 
% joint coordinates for a sequence of points along a trajectory.  Q(i,j) is
% the j'th joint parameter for the i'th trajectory point.  In this case
% FKINE(ROBOT, Q) returns 3D matrix where the last subscript is the index
% along the path.
%
% The robot's base or tool transform, if present, are incorporated into the
% result.
%
% See also: LINK, ROBOT.

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

function t = fkine(robot, q)
	%
	% evaluate fkine for each point on a trajectory of 
	% theta_i or q_i data
	%

	n = robot.n;

	L = robot.link;
	if length(q) == n,
		t = robot.base;
		for i=1:n,
			t = t * L{i}(q(i));
		end
		t = t * robot.tool;
	else
		if numcols(q) ~= n,
			error('bad data')
		end
		t = zeros(4,4,0);
		for qv=q',		% for each trajectory point
			tt = robot.base;
			for i=1:n,
				tt = tt * L{i}(qv(i));
			end
			t = cat(3, t, tt * robot.tool);
		end
	end
