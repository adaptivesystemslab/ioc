%JTRAJ Compute a joint space trajectory between two points
%
%	[Q QD QDD] = JTRAJ(Q0, Q1, N)
%	[Q QD QDD] = JTRAJ(Q0, Q1, N, QD0, QD1)
%	[Q QD QDD] = JTRAJ(Q0, Q1, T)
%	[Q QD QDD] = JTRAJ(Q0, Q1, T, QD0, QD1)
%
% Returns a joint space trajectory Q from state Q0 to Q1.  The number
% of points is N or the length of the given time vector T.  A 7th
% order polynomial is used with default zero boundary conditions for
% velocity and acceleration.  Non-zero boundary velocities can be
% optionally specified as QD0 and QD1.
%
% The function can optionally return a velocity and acceleration
% trajectories as QD and QDD.
%
% Each trajectory is an mxn matrix, with one row per time step, and
% one column per joint parameter.
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

function [qt,qdt,qddt] = jtraj(q0, q1, tv, qd0, qd1)
	if length(tv) > 1,
		tscal = max(tv);
		t = tv(:)/tscal;
	else
		tscal = 1;
		t = [0:(tv-1)]'/(tv-1);	% normalized time from 0 -> 1
	end

	q0 = q0(:);
	q1 = q1(:);

	if nargin == 3,
		qd0 = zeros(size(q0));
		qd1 = qd0;
    elseif nargin == 5,
        qd0 = qd0(:);
        qd1 = qd1(:);
    else
        error('incorrect number of arguments')
	end

	% compute the polynomial coefficients
	A = 6*(q1 - q0) - 3*(qd1+qd0)*tscal;
	B = -15*(q1 - q0) + (8*qd0 + 7*qd1)*tscal;
	C = 10*(q1 - q0) - (6*qd0 + 4*qd1)*tscal;
	E = qd0*tscal; % as the t vector has been normalized
	F = q0;

	tt = [t.^5 t.^4 t.^3 t.^2 t ones(size(t))];
	c = [A B C zeros(size(A)) E F]';
	
	qt = tt*c;

	% compute optional velocity
	if nargout >= 2,
		c = [ zeros(size(A)) 5*A 4*B 3*C  zeros(size(A)) E ]';
		qdt = tt*c/tscal;
	end

	% compute optional acceleration
	if nargout == 3,
		c = [ zeros(size(A))  zeros(size(A)) 20*A 12*B 6*C  zeros(size(A))]';
		qddt = tt*c/tscal^2;
	end
