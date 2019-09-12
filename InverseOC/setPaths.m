% set up the environmental pathing
addpath(genpath('../InverseOC/'));

addpath(genpath('../Utils/')); 
addpath(genpath('../Logic/'));
addpath(genpath('../Libraries/Robotics_Corke/'));
% addpath(genpath('../Libraries/rvctools'));
addpath(genpath('../Libraries/rl'));

% addpath(genpath('../../kalmanfilter/General_FKEKF/DynamicsModelMatlab/MatlabWrapper'));
% addpath(genpath('../../kalmanfilter/ik_framework/instance_iit'));
% addpath(genpath('../../kalmanfilter/ik_framework/instance_jumping'));
% addpath(genpath('../../kalmanfilter/ik_framework/common'));

%% Add Peter Corke's toolbox to library path if required
% if ~exist('rtbdemo', 'file')
%     startup_rvc;
% end