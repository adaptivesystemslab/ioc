function q = angleAxesNotation(rotMtxData)
    if size(rotMtxData, 2) == 1
        theta = rotMtxData;
        w = repmat([0 0 1], length(theta), 1);
        
        q = [theta w];
    else
        theta_prev = 0;
        q = zeros(size(rotMtxData, 1), 4);
        for i = 1:size(rotMtxData, 1)
            currRow = rotMtxData(i, :);
            R = reshape(currRow, [3 3])';

            arg = (trace(R)-1)/2;
            if arg < -1 % numerical issues may push the arg to be slightly less than -1, which will make the acos a complex val
                arg = -1;
            elseif arg > 1
                arg = 1; % same deal here
            end
            
            theta = acos(arg);
            vec = [R(3, 2) - R(2, 3); 
                R(1, 3) - R(3, 1);
                R(2, 1) - R(1, 2)]';
            
            w = (1/(2*sin(theta)))*vec;
            if abs(max(w)) > 1e2 % denom very small number. want to prevent spikes or div0 errors
                % if we get spike, use the offset from the previous
                w = (1/(2*sin(theta_prev)))*vec;
            else
                % there was no spike. we can update the prev theta
                theta_prev = theta;
            end
            
            q(i, :) = [theta w];
        end
        
%         figure; 
%         subplot(211); plot(q(:, 1));
%         subplot(212); plot(q(:, 2:4));
    end
end

function [x,y,z] = decompose_rotationXYZ(R)
	x = atan2(R(3,2), R(3,3));
	y = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
	z = atan2(R(2,1), R(1,1));
end