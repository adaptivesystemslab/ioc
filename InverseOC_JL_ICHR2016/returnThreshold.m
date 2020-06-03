function constThres = returnThreshold(param)
    % tuned to 3 cfs, at 20 slide
    
switch param.win_length
%     case 31
%         constThres = 5.7e-5; % will not work well since it's all constrained
%         
%     case 41
%         constThres = 0; % same here
        
%     case 51
%         constThres = 2.5e-4; % some near misses because of the sim isn't right on a knot point

    case 61
        constThres = 5.0e-3; % 0.005
        
    case 81
        constThres = 6.5e-3;
        
    case 101
        constThres = 4.5e-3;
        
    case 121
        constThres = 3.5e-3;
        
    case 141
        constThres = 3.0e-3;
        
    case 161
        constThres = 2.5e-3;
        
    otherwise
        constThres = 10e-3;    % check1threshold 0.1e-3 1e-3 (1 DOF), 10e-3 (3 DOF)
        
end

% switch param.knots
%     case 0
% end