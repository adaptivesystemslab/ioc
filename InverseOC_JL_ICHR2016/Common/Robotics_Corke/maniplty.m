%MANIPLTY Manipulability measure
%
%	M = MANIPLTY(ROBOT, Q)
%	M = MANIPLTY(ROBOT, Q, WHICH)
%
% Computes the manipulability index for the manipulator at the given pose.
%
% For an n-axis manipulator Q may be an n-element vector, or an m x n
% joint space trajectory.
%
% If Q is a vector MANIPLTY returns a scalar manipulability index.
% If Q is a matrix MANIPLTY returns a column vector of  manipulability 
% indices for each pose specified by Q.
%
% The argument WHICH can be either 'yoshikawa' (default) or 'asada' and
% selects one of two manipulability measures.
% Yoshikawa's manipulability measure gives an indication of how far 
% the manipulator is from singularities and thus able to move and 
% exert forces uniformly in all directions.
%
% Asada's manipulability measure is based on the manipulator's
% Cartesian inertia matrix.  An n-dimensional inertia ellipsoid
% 	X' M(q) X = 1
% gives an indication of how well the manipulator can accelerate
% in each of the Cartesian directions.  The scalar measure computed
% here is the ratio of the smallest/largest ellipsoid axis.  Ideally
% the ellipsoid would be spherical, giving a ratio of 1, but in
% practice will be less than 1.
%
% See also: INERTIA, JACOB0.

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

function w = maniplty(robot, q, which)
	n = robot.n;

	if nargin == 2,
		which = 'yoshikawa';
	end

	if length(q) == robot.n,
		q = q(:)';
	end

	w = [];
	switch which,
	case {'yoshikawa', 'yoshi', 'y'}
		for Q = q',
			w = [w; yoshi(robot, Q)];
		end
	case {'asada', 'a'}
		for Q = q',
			w = [w; asada(robot, Q)];
		end
	end

function m = yoshi(robot, q)
	J = jacob0(robot, q);
	m = sqrt(det(J * J'));

function m = asada(robot, q)
	J = jacob0(robot, q);
	Ji = inv(J);
	M = inertia(robot, q);
	Mx = Ji' * M * Ji;
	e = eig(Mx);
	m = min(e) / max(e);
