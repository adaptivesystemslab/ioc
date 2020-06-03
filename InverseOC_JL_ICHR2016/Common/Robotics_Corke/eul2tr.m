%EUL2TR Convert Euler angles to homogeneous transformation
%
% 	TR = EUL2TR([PHI THETA PSI])
% 	TR = EUL2TR(PHI, THETA, PSI)
%
% Returns a homogeneous tranformation for the specified Euler angles.  These 
% correspond to rotations about the Z, Y, Z axes respectively.
%
% See also: TR2EUL, RPY2TR

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

function T = eul2tr(phi, theta, psi)
	if (nargin == 1),
		if numcols(phi) ~= 3,
			error('bad arguments')
		end
		theta = phi(:,2);
		psi = phi(:,3);
		phi = phi(:,1);
	end

	if numrows(phi) == 1,
                r = rotz(phi) * roty(theta) * rotz(psi);
		T = r2t(r);
	else
		for i=1:numrows(phi),
			r = rotz(phi) * roty(theta) * rotz(psi);
			T(:,:,1) = r2t(r);
		end
	end
