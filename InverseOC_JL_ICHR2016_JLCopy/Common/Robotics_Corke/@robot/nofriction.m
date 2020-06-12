%NOFRICTION Return robot object with zero link friction 
%
%	ROBOT = NOFRICTION(ROBOT)
%
% Return the robot object but with non-linear friction coefficients set to
% zero.  
%
% 	ROBOT = NOFRICTION(ROBOT, 'all')
%
% Return the robot object but with all friction coefficients set to zero.  
%
% Non-linear (Coulomb) friction can cause numerical problems when integrating
% the equations of motion (FDYN).
%
% The resulting robot object has its name string modified by prepending 'NF/'.
%
% See also: LINK/NOFRICTION

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

function  r2 = nofriction(r, varargin)

	r2 = robot(r);

	for i=1:r2.n,
		l2{i} = nofriction(r.link{i}, varargin{:});
	end

	r2.link = l2;
	r2.name = ['NF/' r.name];
