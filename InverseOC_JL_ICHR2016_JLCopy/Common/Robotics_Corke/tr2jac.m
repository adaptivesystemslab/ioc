%TR2JAC Compute a Jacobian to map differentials between frames
%
%	J = TR2JAC(T)
%
% Returns a 6x6 Jacobian matrix to map differentials (joint velocity) between 
% frames related by the homogeneous transform T.
%
% See also: TR2DIFF, DIFF2TR.

% Copyright (C) 1993-2008, by Peter I. Corke
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

function J = tr2jac(t)
	J = [	t(1:3,1)'	cross(t(1:3,4),t(1:3,1))'
		t(1:3,2)'	cross(t(1:3,4),t(1:3,2))'
		t(1:3,3)'	cross(t(1:3,4),t(1:3,3))'
		zeros(3,3)	t(1:3,1:3)'		];
		
