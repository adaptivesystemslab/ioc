%TRINTERP Interpolate homogeneous transformations
%
%	TR = TRINTERP(T0, T1, R)
%
% Returns a homogeneous transform interpolation between T0 and T1 as
% R varies from 0 to 1.  Rotation is interpolated using quaternion
% spherical linear interpolation.
%
% See also: CTRAJ, QUATERNION

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

function t = trinterp(T0, T1, r)

	q0 = quaternion(T0);
	q1 = quaternion(T1);

	p0 = transl(T0);
	p1 = transl(T1);

	qr = qinterp(q0, q1, r);
	pr = p0*(1-r) + r*p1;

	t = [qr.r pr; 0 0 0 1];
