%FRICTION Compute friction torque on the LINK object
%
%	TAU = FRICTION(LINK, QD)
%
% Return the friction torque on the link moving at speed QD.  Depending
% on fields in the LINK object viscous and/or Coulomb friction
% are computed.
%
% SEE ALSO: ROBOT/FRICTION

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

function  tau = friction(l, qd)
	tau = 0.0;

	qd = qd(:);
	tau = l.B * qd;

	tau = tau + (qd > 0) * l.Tc(1) + (qd < 0) * l.Tc(2);
