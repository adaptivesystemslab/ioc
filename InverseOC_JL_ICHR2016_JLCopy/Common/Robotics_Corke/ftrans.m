%FTRANS Transform force/moment
%
%	FT = FTRANS(T, F)
%
% Transforms a force/moment F in the base frame to FT in the frame T.
% F and FT are 6-vectors of the form [Fx Fy Fz Mx My Mz]
%
% SEE ALSO: DIFF2TR

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

function Ft = ftrans(T, F)

	f = F(1:3); m = F(4:6);
	k = cross(f, transl(T) ) + m;

	mt = rot(T)' * k;
	ft = rot(T)' * F(1:3);

	Ft = [ft; mt];
