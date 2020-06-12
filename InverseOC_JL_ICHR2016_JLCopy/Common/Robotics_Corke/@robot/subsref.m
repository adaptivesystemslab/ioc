%SUBSREF Reference methods on a ROBOT object
%
%	ROBOT.n			return number of links
%	ROBOT.link		return cell array of link objects
%	ROBOT.gravity 		return gravity vector
%	ROBOT.base 		return homog xform of robot base
%	ROBOT.tool 		return homog xform of robot tool
%	ROBOT.qlim 		return joint limit matrix
%	ROBOT.offset 		return joint offset vector
%	ROBOT.mdh		return MDH convention boolean (0=DH, 1=MDH)
%
%	ROBOT.islimit 		return joint limit boolean vector
%
%	ROBOT.name		return name of robot
%	ROBOT.manuf		return who built it
%	ROBOT.comment		return general comment
%	ROBOT.config		return joint configuration string
%
%	ROBOT.plotopt 		return options for plot(robot)
%	ROBOT.lineopt 		return line drawing option string for links
%	ROBOT.shadowopt 	return line drawing option string for shadow
%	ROBOT.handle		return graphics handles in object
%	ROBOT.q 		return joint angles for plot(robot)
%
%	ROBOT.dh		return legacy DH matrix
%	ROBOT.dyn		return legacy DYN matrix

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

function v = subsref(r, s)

	if s(1).type  ~= '.'
		%error('only .field supported')
	end

	% NOTE WELL:  the following code can't use getfield() since
	% getfield()  uses this, and Matlab will crash!!

	el = char(s(1).subs);
	switch el,
	%%%%%%%%% retrieve robot parameters
	case 'n',
		v = r.n;
	case 'gravity'
		v = r.gravity;
	case 'tool'
		v = r.tool;
	case 'base'
		v = r.base;
	case 'mdh',
		v = r.mdh;

	case {'link', 'links'},
		if length(s) == 1,
			v = r.link;
		elseif s(2).type == '{}'
			j = s(2).subs;
			j = j{1};
			if (j < 1) | (j > r.n)
				error('link index out of bounds')
			end
			v = r.link{j};
		end
	case 'offset',
		L = r.link;
		v = [];
		for i=1:r.n,
			v = [v; L{i}.offset];
		end
	case 'qlim',
		L = r.link;
		v = [];
		for i=1:r.n,
			v = [v; L{i}.qlim];
		end

	%%%%%%%%% descriptive strings
	case 'name',
		v = r.name;
	case 'manuf',
		v = r.manuf;
	case 'comment',
		v = r.comment;
	case 'config',
		L = r.link;
		v = '';
		for i=1:r.n,
			v(i) = L{i}.RP;
		end

	%%%%%%%%% joint limit test
	case 'islimit',
		L = r.link;
		if s(2).type  ~= '()'
			error('expecting argument for islimit method');
		end
		q = s(2).subs{1};
		if length(q) ~= r.n,
			error('argument for islimit method is wrong length');
		end
		v = [];
		for i=1:r.n,
			v = [v; L{i}.islimit(q(i))];
		end
	%%%%%%%%% legacy DH/DYN support
	case 'dh',
		v = [];
		L = r.link;
		for i=1:r.n,
			v = [v; L{i}.dh];
		end
	case 'dyn'
		v = [];
		L = r.link;
		for i=1:r.n,
			v = [v; L{i}.dyn];
		end

	%%%%%%%%% graphics support
	case 'q',
		v = r.q;
	case 'plotopt',
		v = r.plotopt;
	case 'lineopt'
		v = r.lineopt;
	case 'shadowopt'
		v = r.shadowopt;
	case {'show', 'handle'}
		v = r.handle';
	otherwise, error('Unknown method')
	end
