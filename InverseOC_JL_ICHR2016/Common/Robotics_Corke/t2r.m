%TR2ROT Return rotational submatrix of a homogeneous transformation
%
%	R = TR2ROT(T)
%
% Return R the 3x3 orthonormal rotation matrix from the homogeneous 
% transformation T.
%
% SEE ALSO: ROT2TR

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

function R = tr2rot(T)

	if ~ishomog(T)
		error('input must be a homogeneous transform');
	end

	R = T(1:3,1:3);
