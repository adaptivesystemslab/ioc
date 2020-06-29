%DIFF2TR Convert a differential to a homogeneous transform
%
% 	TR = DIFF2TR(D)
%
% Returns a homogeneous transform representing differential translation 
% and rotation.  The matrix contains a skew symmetric rotation submatrix.
%
% See also: TR2DIFF.

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

function delta = diff2tr(d)
	delta =[	0	-d(6)	d(5)	d(1)
			d(6)	0	-d(4)	d(2)
			-d(5)	d(4)	0	d(3)
			0	0	0	0	];
