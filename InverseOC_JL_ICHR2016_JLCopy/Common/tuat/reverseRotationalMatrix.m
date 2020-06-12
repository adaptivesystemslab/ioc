function q = reverseRotationalMatrix(rotMtxData)
    if size(rotMtxData, 2) == 1
        q = rotMtxData;
    else
        q = zeros(size(rotMtxData, 1), 3);
        for i = 1:size(rotMtxData, 1)
            currRow = rotMtxData(i, :);
            currRot = reshape(currRow, [3 3])';

            % assume XYZ euler
            [x2,y2,z2] = decompose_rotationXYZ(currRot);
            q(i, :) = [x2 y2 z2];
            
%             eul = tr2eul(currRot);   % assuming ZYZs
% %           eul = tr2rpy(currRot);   % assuming rpy
%             q(i, :) = eul;
        end
    end
end

function [x,y,z] = decompose_rotationXYZ(R)
	x = atan2(R(3,2), R(3,3));
	y = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
	z = atan2(R(2,1), R(1,1));
end