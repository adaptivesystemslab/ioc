function [state, control] = setDataLength(state, control, winSize, i)
    lowerBound = i-winSize;
    if lowerBound < 1
        lowerBound = 1;
    end
    indsToUse = lowerBound:i;

    state = state(indsToUse, :);
    control = control(indsToUse, :);
end