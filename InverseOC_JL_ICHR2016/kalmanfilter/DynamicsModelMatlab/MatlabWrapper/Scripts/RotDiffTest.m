%Testing Rotation Differences

a = rand(3,1)*10;
da = rand(3,1)*0.1;

a = [3.1416;   -1.5707;    1.5708];
b = [-0.0000;   -1.5702;   -1.5708];


yaw = a(1);
pitch = a(2);
roll = a(3);

R1 = angle2dcm(yaw,pitch,roll)';
R1z = rotz(yaw);
R1y = roty(pitch);
R1x = rotx(roll);

R2 = angle2dcm(b(1),b(2),b(3))';
R2z = rotz(b(1));
R2y = roty(b(2));
R2x = rotx(b(3));

meanR = SO3.s_exp(SO3.s_hat((SO3.s_hatinv(SO3.s_log(R1))+SO3.s_hatinv(SO3.s_log(R2)))/2));

[y1, p1, r1] = dcm2angle(R1*meanR');
[y2, p2, r2] = dcm2angle(R2*meanR');



%% Testing Lie Group Stuff instead

theta = [rand*10 0 0]';
theta_eps = theta;
theta_eps(1) = theta_eps(1)+1e-6;

genX = [0 0 0; 0 0 -1; 0 1 0];

R01 = SO3.s_exp(SO3.s_hat(rand(3,1)*10));
R12 = SO3.s_exp(SO3.s_hat(rand(3,1)*10));
Rtheta =  SO3.s_exp(SO3.s_hat(theta));
Rtheta_eps =  SO3.s_exp(SO3.s_hat(theta_eps));
R = R01*Rtheta*R12;
Reps = R01*Rtheta_eps*R12;

mes = R*[1 2 3]';
mes_eps = Reps*[1 2 3]';

dmes = (mes_eps-mes)/1e-6;


(SO3.s_hatinv(SO3.s_log(Reps)) - SO3.s_hatinv(SO3.s_log(R)))/1e-6


nw = norm(theta);
Dexp = sin(nw)/nw*eye(3)+ (1-cos(nw))/(nw^2)*skew(theta) +...
    (1-sin(nw)/nw)/(nw^2)*(theta*theta');


%% Test Lie Group Wrapping

pitch = 0:0.01:3*pi;

theta = zeros(numel(pitch),3);
for i=1:numel(pitch)
   r = pitch(i);
   p = 0;
   y = pitch(i);
   R = angle2dcm(y,p,r);
   theta(i,:) = SO3.s_hatinv(SO3.s_log(R));
end

