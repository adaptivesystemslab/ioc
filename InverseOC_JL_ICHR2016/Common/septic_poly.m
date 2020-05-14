function p = septic_poly(q0, dq0, ddq0, dddq0, qf, dqf, ddqf, dddqf, t)
    % 7th order line generator
    
    t0 = t(1);
    tf = t(end);
    
    Q = [1 t0   t0^2   t0^3     t0^4     t0^5      t0^6      t0^7;
         0 1  2*t0   3*t0^2   4*t0^3   5*t0^4    6*t0^5    7*t0^6;
         0 0  2      6*t0    12*t0^2  20*t0^3   30*t0^4   42*t0^5;
         0 0  0      6       24*t0    60*t0^2  120*t0^3  210*t0^4;
         1 tf   tf^2   tf^3     tf^4     tf^5      tf^6      tf^7;
         0 1  2*tf   3*tf^2   4*tf^3   5*tf^4    6*tf^5    7*tf^6;
         0 0  2      6*tf    12*tf^2  20*tf^3   30*tf^4   42*tf^5;
         0 0  0      6       24*tf    60*tf^2  120*tf^3  210*tf^4];
     
    b = [q0,dq0,ddq0,dddq0,qf,dqf,ddqf,dddqf]';
    p = Q\b;
%      a = inv(Q)*[q0,dq0,ddq0,qf,dqf,ddqf]';

%      X = p(1) + p(2)*t + p(3)*t.^2 + p(4)*t.^3   + p(5)*t.^4    + p(6)*t.^5     + p(7)*t.^6     + p(8)*t.^7;
%      dX =       p(2) + 2*p(3)*t  + 3*p(4)*t.^2 + 4*p(5)*t.^3  + 5*p(6)*t.^4   + 6*p(7)*t.^5   + 7*p(8)*t.^6;
%      ddX =             2*p(3)    + 6*p(4)*t   + 12*p(5)*t.^2 + 20*p(6)*t.^3  + 30*p(7)*t.^4  + 42*p(8)*t.^5;
%      dddX =                        6*p(4)     + 24*p(5)*t    + 60*p(6)*t.^2 + 120*p(7)*t.^3 + 210*p(8)*t.^4;
%      ddddX =                                    24*p(5)      +120*p(6)*t    + 360*p(7)*t.^2 + 840*p(8)*t.^3;
     
     p = flipud(p); % flip the a matrix so it matches polyfit
end