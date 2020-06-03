%TRNORM Normalize a homogeneous transformation.
%
%	TN = TRNORM(T) 
%
% Returns a normalized homogeneous tranformation matrix in which the rotation
% submatrix is a proper orthogonal matrix.
% The O and V vectors are normalized and the normal vector is formed from
% O x A.
%
% Finite word length arithmetic can cause transforms to become `unnormalized'.
%
% See also: OA2TR

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

function r = trnorm(t)
	n = cross(t(1:3,2), t(1:3,3));	% N = O x A
    o = cross(t(1:3,3), n);         % O = A x N
	r = [unit(n) unit(t(1:3,2)) unit(t(1:3,3)) t(1:3,4); 0 0 0 1];
