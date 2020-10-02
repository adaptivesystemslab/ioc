%ACCEL Compute manipulator forward dynamics
%
%	QDD = ACCEL(ROBOT, Q, QD, TORQUE)
%	QDD = ACCEL(ROBOT, [Q QD TORQUE])
%
% Returns a vector of joint accelerations that result from applying the 
% actuator TORQUE to the manipulator ROBOT in state Q and QD.
%
% Uses the method 1 of Walker and Orin to compute the forward dynamics.
% This form is useful for simulation of manipulator dynamics, in
% conjunction with a numerical integration function.
%
% See also: RNE, ROBOT, ODE45.

% MOD HISTORY
% 4/99 add object support
% 1/02 copy rne code from inertia.m to here for speed
% $Log: not supported by cvs2svn $
% $Revision: 1.3 $

% Copyright (C) 1999-2008, by Peter I. Corke
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


function qdd = accel(robot, Q, qd, torque)
	n = robot.n;

	if nargin == 2,
	        q = Q(1:n);
		qd = Q(n+1:2*n);
		torque = Q(2*n+1:3*n);
	else
		q = Q;
		if length(q) == robot.n,
			q = q(:);
			qd = qd(:);
		end
	end

	% compute current manipulator inertia
	%   torques resulting from unit acceleration of each joint with
	%   no gravity.
	M = rne(robot, ones(n,1)*q', zeros(n,n), eye(n), [0;0;0]);

	% compute gravity and coriolis torque
	%    torques resulting from zero acceleration at given velocity &
	%    with gravity acting.
	tau = rne(robot, q', qd', zeros(1,n));	

	qdd = inv(M) * (torque(:) - tau');

