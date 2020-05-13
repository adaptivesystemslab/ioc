%CTRAJ Compute a Cartesian trajectory between two points
%
% 	TC = CTRAJ(T0, T1, N)
%	TC = CTRAJ(T0, T1, R)
%
% Returns a Cartesian trajectory TC from point T0 to T1.  The number
% of points is N or the length of the given path distance vector R.
%
% In the first case the points are equally spaced between T0 and T1.
% In the second case R gives the distance along the path, and the 
% elements of R must be in the range [0 1].
%
% Each trajectory is a 4x4xn matrix, with the last subscript being the
% point index.
%
% SEE ALSO: TRINTERP, QINTERP, TRANSL.
%

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

function tt = ctraj(t0, t1, n)
	if length(n) == 1,
		i = 1:n;
		r = (i-1)/(n-1);
	else
		r = n(:)';
		n = length(r);
	end

	if any(r> 1) | any(r<0),
		error('path position values (R) must 0<=R<=1)')
	end
	tt = zeros(4,4,0);

	for R=r,
		tt = cat(3, tt, trinterp(t0, t1, R));
	end

	
