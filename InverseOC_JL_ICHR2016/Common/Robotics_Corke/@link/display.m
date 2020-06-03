%DISPLAY Display the value of a LINK object
%
% This method is invoked to display a link object by the Matlab interpreter,
% and gives a terse single line description of link kinematics.
%
% If invoked with a second argument (value ignored) it will display a long
% form description which includes all defined inertial parameters.

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

function display(l, full)

	disp(' ');
	disp([inputname(1), ' = '])
	disp(' ');
	disp(['  ' char(l)])
	disp(' ');

	if nargin > 1,
		if ~isempty(l.m)
			fprintf('  m    = %f\n', l.m)
		end
		if ~isempty(l.r)
			fprintf('  r    = %f %f %f\n', l.r);
		end
		if ~isempty(l.I)
			fprintf('  I    = | %f %f %f |\n', l.I(1,:));
			fprintf('         | %f %f %f |\n', l.I(2,:));
			fprintf('         | %f %f %f |\n', l.I(3,:));
		end
		if ~isempty(l.Jm)
			fprintf('  Jm   = %f\n', l.Jm);
		end
		if ~isempty(l.B)
			fprintf('  B    = %f\n', l.B);
		end
		if ~isempty(l.Tc)
			fprintf('  Tc   = %f(+) %f(-)\n', l.Tc(1), l.Tc(2));
		end
		if ~isempty(l.G)
			fprintf('  G    = %f\n', l.G);
		end
		if ~isempty(l.qlim)
			fprintf('  qlim = %f to %f\n', l.qlim(1), l.qlim(2));
		end
	end
