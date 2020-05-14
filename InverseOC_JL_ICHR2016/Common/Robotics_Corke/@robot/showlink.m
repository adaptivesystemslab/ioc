%SHOWLINK show all link parameters of ROBOT object
%
%	SHOWLINK(robot)
%
% Shows details of all link parameters for this robot object, including
% inertial parameters.
%
% See also: LINK/SHOWLINK, LINK.

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

function showlink(r)

	l = r.link;
	for i=1:r.n,
		fprintf('Link %d------------------------\n', i);
		showlink(l{i});
	end
