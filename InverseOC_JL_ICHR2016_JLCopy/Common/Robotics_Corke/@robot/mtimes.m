%MTIMES Multiply robot objects
%
% Robot objects can be multiplied r1*r2 which is mechanically equivalent
% to concatenating the two robots, or mounting robot r2 on the end of robot r1.

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

function r2 = mtimes(r, l)

	if ~isa(r, 'robot')
		error('left arg must be a robot')
	end
	if isa(l, 'robot')
		r2 = robot(r);
		r2.link = [r2.link l.link];
		r2.n = length(r2.link);
	elseif isa(l, 'link')
	end
