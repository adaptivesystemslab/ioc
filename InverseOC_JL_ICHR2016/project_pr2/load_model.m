% load model
addpath(genpath('DynamicsModelMatlab'));

mdl = rlCModel('robot.rlmdl.xml');

mdl.bodies(2).m = 25.799322;
mdl.bodies(3).m = 2.74988;
mdl.bodies(4).m = 6.01769;
mdl.bodies(5).m = 1.90327;
mdl.bodies(6).m = 2.57968;
mdl.bodies(7).m = 0.61402;
mdl.bodies(8).m = 0.58007;
mdl.bodies(9).m = 0.17126;
mdl.bodies(10).m = 0.04419;
mdl.bodies(11).m = 0.17389;
mdl.bodies(12).m = 0.04419;

% load data
process_data