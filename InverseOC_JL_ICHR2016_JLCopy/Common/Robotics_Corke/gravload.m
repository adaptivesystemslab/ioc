%GRAVLOAD Compute the gravity loading on manipulator joints
%
%	TAUG = GRAVLOAD(ROBOT, Q)
%	TAUG = GRAVLOAD(ROBOT, Q, GRAV)
%
% Compute the joint gravity loading for the manipulator ROBOT in the
% configuration Q.
%
% If Q is a row vector, the result is a row vector of joint torques.
% If Q is a matrix, each row is interpretted as a joint state vector, and
% the result is a matrix each row being the corresponding joint torques.
%
% Gravity vector can be given explicitly using the GRAV argument, otherwise
% it defaults to the value of the ROBOT object.
%
% See also: ROBOT, RNE, ITORQUE, CORIOLIS.

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

function tg = gravload(robot, q, grav)
	if numcols(q) ~= robot.n
		error('Insufficient columns in q')
	end
	if nargin == 2,
		tg = rne(robot, q, zeros(size(q)), zeros(size(q)));
	elseif nargin == 3,
		tg = rne(robot, q, zeros(size(q)), zeros(size(q)), grav);
	end
	
