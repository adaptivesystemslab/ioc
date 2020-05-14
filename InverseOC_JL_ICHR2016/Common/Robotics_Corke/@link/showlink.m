%SHOWLINK show all parameters of LINK object
%
%	SHOWLINK(link)
%
% Display all parameters, including all defined inertial parameters, for the
% link object.
%
% See also: ROBOT/SHOWLINK, LINK.

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

function showlink(l)

	llab = 6;
	for n =fieldnames(l)'
		v = getfield(l, char(n));
		name = char(n);
		spaces = char(' '*ones(1,llab-length(name)));
		if min(size(v)) == 1,
			val = num2str(v(:)');
		else
			val = num2str(v);
		end
		label = [name spaces ' = '];
		if numrows(val) > 1,
			pad = {label; char(' '*ones(numrows(val)-1,1))};
		else
			pad = label;
		end
		disp([char(pad) val]);
	end
