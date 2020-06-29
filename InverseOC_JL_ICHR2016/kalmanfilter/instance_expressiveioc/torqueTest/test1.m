% testing
% q = [pi/2 pi/3 pi/4 pi/2 pi/6];
q = [0 0 0 0 0];

twolink.plot(q);

mdl.position = q;
mdl.forwardPosition();
vis.update();