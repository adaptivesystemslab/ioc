%R2T Return a homogeneous transformation from a rotation
%
%	T = R2T(R)
%
% Return T the 4x4 the homogeneous transformation from the 3x3 rotation R
%
% SEE ALSO: T2R

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

function T = r2t(R)

	T = [R(1:3,1:3) [0;0;0]; 0 0 0 1];
