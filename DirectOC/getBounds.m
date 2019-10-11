function [state, control] = getBounds(motion)
% Get min and max value of each variable in state and control vectors
    state.low = [min(motion.q), min(motion.dq)]';
    state.upp = [max(motion.q), max(motion.dq)]';

    control.low = min(motion.tau)';
    control.upp = max(motion.tau)';
end

