%FDYN Integrate forward dynamics
%
%	[T Q QD] = FDYN(ROBOT, T0, T1)
%	[T Q QD] = FDYN(ROBOT, T0, T1, TORQFUN)
%	[T Q QD] = FDYN(ROBOT, T0, T1, TORQFUN, Q0, QD0)
%	[T Q QD] = FDYN(ROBOT, T0, T1, TORQFUN, Q0, QD0, ARG1, ARG2, ...)
%
% Integrates the dynamics of manipulator ROBOT dynamics over the time 
% interval T0 to T1 and returns vectors of joint position and velocity.
% ROBOT is a robot object and describes the manipulator dynamics and 
% kinematics, and Q is an n element vector of joint state.
%
% A control torque may be specified by a user specified function
%
% 	TAU = TORQFUN(T, Q, QD, ARG1, ARG2, ...)
%
% where Q and QD are the manipulator joint coordinate and velocity state 
% respectively], and T is the current time. Optional arguments passed to FDYN
% will be passed through to the user function.
%
% If TORQFUN is not specified, or is given as 0,  then zero torque is 
% applied to the manipulator joints.
%
% See also: ACCEL, NOFRICTION, RNE, ROBOT, ODE45.

% Copyright (C) 1993-2008 Peter Corke
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

function [t, q, qd] = fdyn(robot, t0, t1, torqfun, q0, qd0, varargin)

	% check the Matlab version, since ode45 syntax has changed
	v = ver;
	if str2num(v(1).Version)<6,
		%error('fdyn now requires Matlab version >= 6');
	end

	n = robot.n;
	if nargin == 3,
		torqfun = 0;
		x0 = zeros(2*n,1);
	elseif nargin == 4,
		x0 = zeros(2*n, 1);
	elseif nargin >= 6,
		x0 = [q0(:); qd0(:)];
	end
		
	[t,y] = ode45('fdyn2', [t0 t1], x0, [], robot, torqfun, varargin{:});
	q = y(:,1:n);
	qd = y(:,n+1:2*n);

