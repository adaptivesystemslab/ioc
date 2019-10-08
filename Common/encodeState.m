function state = encodeState(position, velocity)
%     lenState = size(position, 2)*2;
%     indPos = 1:2:lenState;
%     indVel = 2:2:lenState;
%     state(:, indPos) = position;
%     state(:, indVel) = velocity;  

    state = [position velocity];
end