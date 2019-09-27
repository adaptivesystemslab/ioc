addpath(genpath('..\Libraries\rl'));

dt = 0.01;
y = [0 0 0 0;
    60 -45 45 90];
dy = [0 0 0 0;
    0 0 0 0];
ddy = [0 0 0 0;
    0 0 0 0];
x = [1
    100];

% generate quintic spline
[q, dq, ddq, t] = quinticSpline(deg2rad(y), dy, ddy, x, dt);

% now use RL to generate tau
tau = zeros(size(q));
rlModel = define4DofArm('rl');
for i = 1:size(q, 1)
    currQ = q(i, :);
    currDq = dq(i, :);
    currDdq = ddq(i, :);
    
    rlModel.updateState(currQ(1:end), currDq(1:end));
    tau(i, :) = rlModel.inverseDynamics(currDdq);
end

figure; 
subplot(221);
plot(t, q); title('q');
hold on
for i = size(x)
    plot(t(x([i i])), [-1 1]);
end

subplot(222);
plot(t, dq); title('dq');
subplot(223);
plot(t, ddq); title('ddq');
subplot(224);
plot(t, tau); title('tau');