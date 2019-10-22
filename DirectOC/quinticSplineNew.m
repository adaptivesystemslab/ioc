function [q, dq, ddq, time] = quinticSplineNew(y, dy, ddy, t, dt)  
%	generates a sequence of quintic splines based on via points. it does so
%	by using the Vandermonde matrix to generate the spline polynomial
%	coefficients, then use polyval to fill in the spline

% 	y(viapoints, dof) - this determines the key poses that the spline needs
%       to pass through. The first dim is the different key poses, while 
%       the second is the different dofs. A similar structure holds for dy 
%       and ddy
%   x(index) - this determines the time index of where the said via point
%       will occur. That is, the last index in x denotes the total length
%       of the spline
%   dt - the delta time of the data

    q=[];
    dq=[];
    ddq=[];
    time=[];
    for ind_y = 1:size(y, 1)-1
        t_curr=t(ind_y):dt: t(ind_y+1);
        seg_q=[];
        seg_dq=[];
        seg_ddq=[];
        for ind_dof = 1:size(y, 2)
            sp_qd0 = quintic_poly(y(ind_y, ind_dof), dy(ind_y, ind_dof), ddy(ind_y, ind_dof), ...
                y(ind_y+1, ind_dof), dy(ind_y+1, ind_dof), ddy(ind_y+1, ind_dof), ...
                t_curr);
            sp_qd1 = polyder(sp_qd0);
            sp_qd2 = polyder(sp_qd1);
            
            seg_q(:,ind_dof)=polyval(sp_qd0, t_curr);
            seg_dq(:,ind_dof)=polyval(sp_qd1, t_curr);
            seg_ddq(:,ind_dof)=polyval(sp_qd2, t_curr);
        end
        q(end+1:end+length(t_curr),:)=seg_q;
        dq(end+1:end+length(t_curr),:)=seg_dq;
        ddq(end+1:end+length(t_curr),:)=seg_ddq;
        time(end+1:end+length(t_curr))=t_curr;
    end
end

function p = quintic_poly(q0, dq0, ddq0, qf, dqf, ddqf, t)
    % 5th order spline, using the Vandermonde method
    t0 = t(1);
    tf = t(end);
    
    Q = [1 t0   t0^2   t0^3     t0^4     t0^5;
         0 1  2*t0   3*t0^2   4*t0^3   5*t0^4;
         0 0  2      6*t0    12*t0^2  20*t0^3;
         1 tf   tf^2   tf^3     tf^4     tf^5;
         0 1  2*tf   3*tf^2   4*tf^3   5*tf^4;
         0 0  2      6*tf    12*tf^2  20*tf^3];
    b = [q0,dq0,ddq0,qf,dqf,ddqf]';
    
    p = Q\b;
    
%      a = inv(Q)*[q0,dq0,ddq0,qf,dqf,ddqf]';

%      X = p(1) + p(2)*t + p(3)*t.^2 + p(4)*t.^3   + p(5)*t.^4    + p(6)*t.^5;
%      dX =       p(2) + 2*p(3)*t  + 3*p(4)*t.^2 + 4*p(5)*t.^3  + 5*p(6)*t.^4;
%      ddX =             2*p(3)    + 6*p(4)*t   + 12*p(5)*t.^2 + 20*p(6)*t.^3;
%      dddX =                        6*p(4)     + 24*p(5)*t    + 60*p(6)*t.^2;
     
     p = flipud(p)'; % flip the a matrix so it matches polyfit
end