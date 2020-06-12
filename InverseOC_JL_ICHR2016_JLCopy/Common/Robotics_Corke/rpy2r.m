%RPY2R Roll/pitch/yaw to rotation matrix
%
% 	R = RPY2R([R P Y])
%	R = RPY2R(R,P,Y)
%
% Returns a 3x3 rotation matrix for the specified roll/pitch/yaw angles.
% These correspond to rotations about the X, Y, Z axes respectively.
%
% NOTE: in previous releases (<8) the angles corresponded to rotations about Z, Y, X.
%
% See also: TR2RPY, EUL2TR

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

function r = rpy2r(roll, pitch, yaw)
	if (nargin == 1),
		if numcols(roll) ~= 3,
			error('bad arguments')
		end
		pitch = roll(:,2);
		yaw = roll(:,3);
		roll = roll(:,1);
	end

	if numrows(roll) == 1,
		r = rotx(roll) * roty(pitch) * rotz(yaw);
	else
		for i=1:numrows(roll),
			r(:,:,i) = rotx(roll(i)) * roty(pitch(i)) * rotz(yaw(i));
		end
	end
