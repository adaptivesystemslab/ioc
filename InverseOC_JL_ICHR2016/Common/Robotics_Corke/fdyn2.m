%FDYN2  private function called by FDYN
%
%	XDD = FDYN2(T, X, FLAG, ROBOT, TORQUEFUN)
%
% Called by FDYN to evaluate the robot velocity and acceleration for
% forward dynamics.  T is the current time, X = [Q QD] is the state vector,
% ROBOT is the object being integrated, and TORQUEFUN is the string name of
% the function to compute joint torques and called as
%
%       TAU = TORQUEFUN(T, X)
%
% if not given zero joint torques are assumed.
%
% The result is XDD = [QD QDD].
%
% See also: FDYN

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

function xd = fdyn2(t, x, flag, robot, torqfun, varargin)

	n = robot.n;

	q = x(1:n);
	qd = x(n+1:2*n);

	% evaluate the torque function if one is given
	if isstr(torqfun)
		tau = feval(torqfun, t, q, qd, varargin{:});
	else
		tau = zeros(n,1);
	end
	
	qdd = accel(robot, x(1:n,1), x(n+1:2*n,1), tau);
	xd = [x(n+1:2*n,1); qdd];
