function [position, velocity] = decodeState(state)
%     lenState = size(state, 2);
%     position = state(:, 1:2:lenState);
%     velocity = state(:, 2:2:lenState);

    lenPos = size(state, 2)/2;
    position = state(:, 1:lenPos);
    velocity = state(:, (lenPos+1):end);
end