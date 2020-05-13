%SUBSREF Reference methods on a QUATERNION object
%
%	QUATERNION.d		return a 4-vector of quaternion elements
%	QUATERNION.s		return the scalar component
%	QUATERNION.v		return the vector component
%	QUATERNION.t		return a 4x4 homogeneous transform
%	QUATERNION.r		return a 3x3 orthonormal rotation matrix

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

function v = subsref(q, s)
	if s(1).type  == '.'

		% NOTE WELL:  the following code can't use getfield() since
		% getfield()  uses this, and Matlab will crash!!

		el = char(s(1).subs);
		switch el,
		case 'd',
			v = double(q);
		case 's',
			v = q.s;
		case 'v',
			v = q.v;
		case 't',
			v = q2tr(q);
		case 'r',
			v = q2tr(q);
			v = v(1:3,1:3);
		end
	else
		error('only .field supported')
	end

%Q2TR	Convert unit-quaternion to homogeneous transform
%
%	T = q2tr(Q)
%
%	Return the rotational homogeneous transform corresponding to the unit
%	quaternion Q.
%
%	See also: TR2Q

%	Copyright (C) 1993 Peter Corke
function t = q2tr(q)

	q = double(q);
	s = q(1);
	x = q(2);
	y = q(3);
	z = q(4);

	r = [	1-2*(y^2+z^2)	2*(x*y-s*z)	2*(x*z+s*y)
		2*(x*y+s*z)	1-2*(x^2+z^2)	2*(y*z-s*x)
		2*(x*z-s*y)	2*(y*z+s*x)	1-2*(x^2+y^2)	];
	t = eye(4,4);
	t(1:3,1:3) = r;
	t(4,4) = 1;
