%TRPLOT Plot a transformation
%
%   TRPLOT(T)
%   TRPLOT(T, name)
%   TRPLOT(T, name, color)
%
% In a set of axes draw a coordinate frame.  The frame can optionally
% be named and have a specified color.

% Copyright (C) 2008, by Peter I. Corke
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

function trplot(T, name, color)

	if nargin == 1,
		fmt = '%c';
	else
		fmt = sprintf('%%c%s', name);
	end
	if nargin < 3,
		color = 'b';
	end

	q = quaternion(T);
	plot(q, transl(T), fmt, color);
