%NOFRICTION Return link object with zero friction 
%
%	LINK = NOFRICTION(LINK)
%	LINK = NOFRICTION(LINK, 'all')
%
% Return the link object with Coulomb or all friction terms set to zero.
%
% See also: ROBOT/NOFRICTION

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

function  l2 = nofriction(l, only)

	l2 = link(l);

	if (nargin == 2) & strcmpi(only(1:3), 'all'),
		l2.B = 0;
	end
	l2.Tc = [0 0];
