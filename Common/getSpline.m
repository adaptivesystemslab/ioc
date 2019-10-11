function [angles, velocities, control, time] = getSpline(model,...
    initialState, finalState, dt)
    
    % This function generates a quantic spline as the initial guess of a
    % trajectory optimization instance

    if ~exist('dt', 'var')
        dt = 0.01;
    end
    
    y = [initialState.jointAngles';
        finalState.jointAngles'];
    
    dy = [initialState.angularVelocities';
          finalState.angularVelocities'];
    
    ddy = zeros(2,length(initialState.jointAngles));
    
    x = [1 100];
    
    [angles, velocities, accelerations, time] = quinticSpline(deg2rad(y),...
        dy, ddy, x, dt);
    
    control = zeros(size(angles));
    
    previousState = model.getState();
    
    for i = 1:size(angles, 1)
        currQ = angles(i, :);
        currDq = velocities(i, :);
        currDdq = accelerations(i, :);
    
        model.updateState(currQ(1:end), currDq(1:end));
        control(i, :) = model.inverseDynamics(currDdq);
    end
  
    model.updateState(previousState(1,:), previousState(2,:));
    angles = rad2deg(angles);
    
%     figure; 
%     subplot(221);
%     plot(time, angles); title('q');
%     hold on
%     for i = size(x)
%         plot(time(x([i i])), [-1 1]);
%     end
% 
%     subplot(222);
%     plot(time, velocities); title('dq');
%     subplot(223);
%     plot(time, accelerations); title('ddq');
%     subplot(224);
%     plot(time, control); title('tau');
end

