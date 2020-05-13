%SUBSASGN assignment methods on a LINK object
%
%	LINK.alpha = alpha	kinematic parameters
%	LINK.A = A
%	LINK.theta = theta
%	LINK.D = D
%	LINK.sigma = sigma	1 if joint is prismatic
%
%	LINK.I = 3x3 inertia matrix about link COG
%	LINK.I = 6x1 inertia vector [Ixx Iyy Izz Ixy Iyz Ixz] about link COG
%	LINK.m = link mass
%	LINK.r = 3vector  link COG wrt link coordinate frame
%
%	LINK.B = link viscous friction (motor referred)
%	LINK.Tc = link Coulomb friction 1 element if symmetric, else 2
%
%	LINK.G = G		gear ratio
%	LINK.Jm = Jm		gear ratio (motor referred)
%
%	LINK.qlim = 2 x 1	joint limit matrix [lower upper]
%	LINK.offset = q0	joint coordinate offset

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

function l = subsasgn(l, s, v)

	if s(1).type  ~= '.'
		error('only .field supported')
	end
	switch s(1).subs,
	case 'alpha',
		l.alpha = v;
	case 'A',
		l.A = v;
	case 'theta',
		l.theta = v;
	case 'D',
		l.D = v;
	case 'offset',
		l.offset = v;
	case 'sigma',
		if ischar(v)
			l.sigma = lower(v) == 'p';
		else
			l.sigma = v;
		end
	case 'G',
		l.G = v;
	case 'I',
		if all(size(v) == [3 3])
			l.I = v;
		elseif length(v) == 3,
			l.I = diag(v);
		elseif length(v) == 6,
			l.I = [	v(1) v(4) v(6)
				v(4) v(2) v(5)
				v(6) v(5) v(3)];
		end
	case 'r',
		l.r = v(:);	% a column vector
	case 'Jm',
		l.Jm = v;
	case 'B',
		l.B = v;
	case 'Tc',
		if length(v) == 1,
			l.Tc = [v -v];
		elseif length(v) == 2,
			l.Tc = v;
		else
			error('Coulomb friction vector can have 1 (symmetric) or 2 (asymmetric) elements only')
		end
	case 'm',
		l.m = v;
	case 'qlim',
		if length(v) ~= 2,
			error('joint limit must have 2 elements');
		end
		l.qlim = v;
	otherwise, error('Unknown method')
	end
