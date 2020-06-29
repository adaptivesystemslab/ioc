function p = quinticPoly(q0, dq0, ddq0, qf, dqf, ddqf, t0, tf)
    % 5th order line generator
    
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