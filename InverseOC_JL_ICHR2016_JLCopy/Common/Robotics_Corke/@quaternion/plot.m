%PLOT Plot a quaternion object 
%
%	PLOT(Q)
%
% Display the quaternion as a rotated coordinate frame.
%
% SEE ALSO: QUATERNION

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

function plot(Q, off, fmt, color)
	%axis([-1 1 -1 1 -1 1])

	if nargin < 2,
		off = [0 0 0];
	end
	if nargin < 3,
		fmt = '%c';
	end
	if nargin < 4,
		color = 'b';
	end
	% create unit vectors
	o = [0 0 0]';
	x1 = Q*[1 0 0]';
	y1 = Q*[0 1 0]';
	z1 = Q*[0 0 1]';

	get(gca, 'Tag')
	if strcmp(get(gca, 'Tag'), 'trplot') == 0,
		fprintf('No tag\n');
		clf
		axes
		set(gca, 'Tag', 'trplot')
		fprintf('set tag\n');
		xlabel( 'X');
		ylabel( 'Y');
		zlabel( 'Z');
	end
	ih = ishold;
	hold on
	plot3([0;x1(1)]+off(1), [0; x1(2)]+off(2), [0; x1(3)]+off(3), color);
	h = text(off(1)+x1(1), off(2)+x1(2), off(3)+x1(3), sprintf(fmt, 'X'));
	set(h, 'Color', color);

	plot3([0;y1(1)]+off(1), [0; y1(2)]+off(2), [0; y1(3)]+off(3), color);
	h = text(off(1)+y1(1), off(2)+y1(2), off(3)+y1(3), sprintf(fmt, 'Y'));
	set(h, 'Color', color);

	plot3([0;z1(1)]+off(1), [0; z1(2)]+off(2), [0; z1(3)]+off(3), color);
	h = text(off(1)+z1(1), off(2)+z1(2), off(3)+z1(3), sprintf(fmt, 'Z'));
	set(h, 'Color', color);
	grid on
	if ~ishold,
		hold off
	end
	axis equal
