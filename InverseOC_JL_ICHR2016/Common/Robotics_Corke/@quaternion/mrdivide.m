%MRDIVIDE Compute quaternion quotient.
%
% Invoked on the / operator, handle two cases:
% q1/q2  	multiply one quaternion by inverse of the second.
% q1/s		result is non-unit quaternion, all elements divided by s

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

function qq = mrdivide(q1, q2)

	if isa(q2, 'quaternion'),
		% qq = q1 / q2
		%    = q1 * qinv(q2)

		qq = q1 * inv(q2);
	elseif isa(q2, 'double'),
		qq = quaternion( double(q1) / q2 );
	end
