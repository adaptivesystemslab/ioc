%ROBOT Robot object constructor
%
%	ROBOT			create a ROBOT object with no links
%	ROBOT(robot)		create a copy of an existing ROBOT object
%	ROBOT(robot, LINK)	replaces links for robot object
%	ROBOT(LINK, ...)	create from a cell array of LINK objects
%	ROBOT(DH, ...)		create from legacy DH matrix
%	ROBOT(DYN, ...)		create from legacy DYN matrix
%
% Optional trailing arguments are:
% 	Name			robot type or name
% 	Manufacturer		who built it
% 	Comment			general comment
%
% If the legacy matrix forms are used the default name is the workspace
% matrix that held the data.
%
% See also: ROBOT/SUBSREF, ROBOT/SUBSASGN, LINK.

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

function r = robot(L, a1, a2, a3)

	if nargin == 0
		r.name = 'noname';
		r.manuf = '';
		r.comment = '';
		r.link = {};
		r.n = 0;
		r.mdh = 0;
		r.gravity = [0; 0; 9.81];
		r.base = eye(4,4);
		r.tool = eye(4,4);
		r.handle = [];	% graphics handles
		r.q = [];	% current joint angles
		r.plotopt = {};
		r.lineopt = {'Color', 'black', 'Linewidth', 4};
		r.shadowopt = {'Color', 'black', 'Linewidth', 1};
		r = class(r, 'robot');
	elseif isa(L, 'robot')
		r = L;
		if nargin == 2,
			r.link = a1;
		end
	else
		% assume arguments are: name, manuf, comment
		if nargin > 1,
			r.name = a1;
		else
			r.name = 'noname';
		end
		if nargin > 2,
			r.manuf = a2;
		else
			r.manuf = '';
		end
		if nargin > 3,
			r.comment = a3;
		else
			r.comment = '';
		end

		if isa(L, 'double')
			% legacy matrix
			dh_dyn = L;
			clear L
			for j=1:numrows(dh_dyn)
				L{j} = link(dh_dyn(j,:));
			end
			% get name of variable
			r.name = inputname(1);
			r.link = L;
		elseif iscell(L) & isa(L{1}, 'link')

			r.link = L;
		else
			error('unknown type passed to robot');
		end
		r.n = length(L);

		% set the robot object mdh status flag
		mdh = [];
		for j = 1:length(L)
			mdh = [mdh L{j}.mdh];
		end
		if all(mdh == 0)
			r.mdh = mdh(1);
		elseif all (mdh == 1)
			r.mdh = mdh(1);
		else
			error('robot has mixed D&H link conventions');
		end

		% fill in default base and gravity direction
		r.gravity = [0; 0; 9.81];
		r.base = eye(4,4);
		r.tool = eye(4,4);
		r.handle = [];
		r.q = [];
		r.plotopt = {};
		r.lineopt = {'Color', 'black', 'Linewidth', 4};
		r.shadowopt = {'Color', 'black', 'Linewidth', 1};

		r = class(r, 'robot');
	end
