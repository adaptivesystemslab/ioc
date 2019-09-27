% Parent function to tune process noise of EKF for jumping data
clear all;

procNoise_lb = 0.01; % this is observation noise, shouldn't be less than this
procNoise_ub = 100; % seems very large, shouldn't be more than this
options = optimoptions('fmincon','Display','off');

tic;
procNoise_optim = fmincon(@tuneEKFParams_allParts,1,[],[],[],[],procNoise_lb,procNoise_ub,[],options)
endTime = toc

