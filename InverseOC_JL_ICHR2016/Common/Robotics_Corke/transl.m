%TRANSL Create translational transform
%
%	TR = TRANSL(X, Y, Z)
%	TR = TRANSL( [X Y Z] )
%
% Returns a homogeneous transformation representing a translation of X, Y
% and Z.
%
%	[X Y Z]' = TRANSL(T)
%
% Returns the translational part of a homogenous transform as a 3-element 
% column vector.
%
%	[X Y Z] = TRANSL(TG)
%
% Returns a  matrix of the X, Y and Z elements extracted from a Cartesian 
% trajectory matrix TG.
%
% See also: CTRAJ.

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

function r = transl(x, y, z)
	if nargin == 1,
		if ishomog(x),
			r = x(1:3,4);
		elseif ndims(x) == 3,
			r = squeeze(x(1:3,4,:))';
		else
			t = x(:);
			r =    [eye(3)			t;
				0	0	0	1];
		end
	elseif nargin == 3,
		t = [x; y; z];
		r =    [eye(3)			t;
			0	0	0	1];
	end
