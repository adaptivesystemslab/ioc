%OA2R Convert O/A vectors to rotation matrix
%
% 	R = OA2R(O, A)
%
% Returns a rotation matrix for the specified orientation and 
% approach vectors formed from 3 vectors such that
% R = [N O A] and N = O x A.  
% The submatrix is guaranteed to be orthonormal so long as O and A are 
% not parallel.
%
% See also: RPY2R, EUL2R, OA2TR

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

function r = oa2r(o, a)
	n = cross(o, a);
    o = cross(a, n);
	r = [unit(n(:)) unit(o(:)) unit(a(:))];
