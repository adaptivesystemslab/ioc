function [q, polyParam] = calc_trajectory_septic(t, y)
    Nh = floor(length(t)/2);
    Nf = length(t);
    
    for j = 1:size(y, 1)
        % only the start and end is passed in, use the 7th order
        % polynomial fit
        [y_start, dy_start, ddy_start, dddy_start, ddddy_start, a_start] = Septic(y(j, 1)  ,0,0,0,y(j, end)  ,0,0,0,t(1:Nh));
%                 [y_end, dy_end, ddy_end, dddy_end, ddddy_end, a_end]             = Septic(y(j, 2)  ,0,0,0,y(j, 1)  ,0,0,0,t(Nh+1:Nf));

%                 q(j, :) = [y_start y_end];
%                 dq(j, :) = [dy_start dy_end];
%                 ddq(j, :) = [ddy_start ddy_end];
%                 dddq(j, :) = [dddy_start dddy_end];
        
        q(j, :) = [y_start fliplr(y_start)];
        polyParam(j, :) = a_start;
    end
end