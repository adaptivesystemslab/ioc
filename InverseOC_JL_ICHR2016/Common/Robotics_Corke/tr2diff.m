%TR2DIFF Convert a transform difference to differential representation
%
%	D = TR2DIFF(T)
%	D = TR2DIFF(T1, T2)
%
% First form converts a homogeneous transform representing an
% infinitessimal motion to a 6-element differential representation.
% Such a homogeneous transform has a rotational submatrix that is,
% approximately, skew symmetric.
%
% Second form returns the 6-element differential motion required to move
% from T1 to T2 in base coordinates.
%
% See also: DIFF2TR.

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

function d = tr2diff(t1, t2)
	if nargin == 1,
		d = [	t1(1:3,4);
			0.5*[t1(3,2)-t1(2,3); t1(1,3)-t1(3,1); t1(2,1)-t1(1,2)]];
	else
		d = [	t2(1:3,4)-t1(1:3,4);
			0.5*(	cross(t1(1:3,1), t2(1:3,1)) + ...
				cross(t1(1:3,2), t2(1:3,2)) + ...
				cross(t1(1:3,3), t2(1:3,3)) ...
			)];
	end

