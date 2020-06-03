%CORIOLIS Compute the manipulator Coriolis matrix
%
% 	C = CORIOLIS(ROBOT, Q, QD)
%
% Returns the n element Coriolis/centripetal torque vector at the specified 
% pose and velocity.
% ROBOT is a robot object and describes the manipulator dynamics and 
% kinematics.
%
% If Q and QD are row vectors, CORIOLIS(ROBOT,Q,QD) is a row vector 
% of joint torques.
% If Q and QD are matrices, each row is interpretted as a joint state 
% vector, and CORIOLIS(ROBOT,Q,QD) is a matrix each row being the 
% corresponding joint %	torques.
%
% See also: ROBOT, RNE, ITORQUE, GRAVLOAD.


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

function c = coriolis(robot, q, qd)

	if nargin == 3,
		c = rne(robot, q, qd, zeros(size(q)), [0;0;0]);
	else
		n = length(q);
		c = [];
		qd = zeros(1,n);
		for i=1:n,
			qd(i) = 1;
			C = coriolis(robot, q, qd);
			qd(i) = 0;
			c(:,i) = C';
		end
	end
