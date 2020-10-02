%CHAR String representation of robot parametesrs

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

function s = char(r)

	% build a configuration string
	rp = [];
	for i = 1:r.n,
		rp = [rp r.link{i}.RP];
	end

	s = sprintf('%s (%d axis, %s)', r.name, r.n, rp);

	if ~isempty(r.manuf)
		s = strcat(s, [' [' r.manuf ']']);
	end
	if ~isempty(r.comment)
		s = strcat(s, [' <' r.comment '>']);
	end
	s = strcat(s, sprintf('\n\t\tgrav = [%.2f %.2f %.2f]\n', r.gravity));
	if getfield(r, 'mdh') == 0,
		s = strcat(s, sprintf('\t\tstandard D&H parameters\n'));
	else
		s = strcat(s, sprintf('\t\tmodified D&H parameters\n'));
	end

	s = strcat(s, sprintf('\n\n  alpha\t\t A\t\t theta\t\t D\t\tR/P\n'));
	for i = 1:r.n,
		s = strcat(s, sprintf('\n%s', char(r.link{i})));
	end
