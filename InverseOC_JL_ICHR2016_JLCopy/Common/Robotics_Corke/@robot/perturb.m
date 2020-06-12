%PERTURB Return robot object with perturbed dynamic parameters
%
%	ROBOT = PERTURB(ROBOT, P)
%
% Return a new robot object in which the dynamic parameters (link mass and
% inertia) have been perturbed.  The perturbation is multiplicative so that
% values are multiplied by random numbers in the interval (1-P) to (1+P).
%
% Useful for investigating the robustness of various model-based control 
% schemes.
%
% The name string of the perturbed robot is prefixed by 'P/'.
%

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

function  r2 = perturb(r, p)

	if nargin == 1,
		p = 0.1;	% 10 percent disturb by default
	end


	for i=1:r.n,
		l2{i} = r.link{i};
		s = (2*rand-1)*p + 1;
		l2{i}.m = l2{i}.m * s;
		s = (2*rand-1)*p + 1;
		l2{i}.I = l2{i}.I * s;
	end

	r2 = robot(r, l2);		% clone the robot
	r2.name = ['P/' r.name];
