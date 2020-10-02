%JACOBN Compute manipulator Jacobian in end-effector frame
%
%	JN = JACOBN(ROBOT, Q)
%
% Returns a Jacobian matrix for the robot ROBOT in pose Q.
%
% The manipulator Jacobian matrix maps differential changes in joint space
% to differential Cartesian motion of the end-effector (end-effector coords).
% 		dX = J dQ
%
% This function uses the technique of
% 	Paul, Shimano, Mayer
% 	Differential Kinematic Control Equations for Simple Manipulators
% 	IEEE SMC 11(6) 1981
% 	pp. 456-460
%
% For an n-axis manipulator the Jacobian is a 6 x n matrix.
%
% See also: JACOB0, DIFF2TR, TR2DIFF

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

function J = jacobn(robot, q)

	n = robot.n;
	L = robot.link;		% get the links

	J = [];
	U = robot.tool;
	for j=n:-1:1,
		if robot.mdh == 0,
			% standard DH convention
			U = L{j}( q(j) ) * U;
		end
		if L{j}.RP == 'R',
			% revolute axis
			d = [	-U(1,1)*U(2,4)+U(2,1)*U(1,4)
				-U(1,2)*U(2,4)+U(2,2)*U(1,4)
				-U(1,3)*U(2,4)+U(2,3)*U(1,4)];
			delta = U(3,1:3)';	% nz oz az
		else
			% prismatic axis
			d = U(3,1:3)';		% nz oz az
			delta = zeros(3,1);	%  0  0  0
		end
		J = [[d; delta] J];

		if robot.mdh ~= 0,
			% modified DH convention
			U = L{j}( q(j) ) * U;
		end
	end
