% set up the environmental pathing
addpath(genpath('./Support/'));

addpath(genpath('../Utils/')); 
addpath(genpath('../Logic/'));
addpath(genpath('../Libraries/rvctools'));
addpath(genpath('../Libraries/rl'));

addpath(genpath('../../kalmanfilter/General_FKEKF/DynamicsModelMatlab/MatlabWrapper'));
addpath(genpath('../../kalmanfilter/ik_framework/instance_iit'));
addpath(genpath('../../kalmanfilter/ik_framework/common'));

% name = getenv('COMPUTERNAME');
% switch name
%     case 'MSI-JLIN'
%        
%     otherwise
%         
% end

%  addpath('../Libraries/cvx');
 
% switch isunix
%     case 0
%         addpath('../Libraries/cvx_w64');
%         
%     case 1
%         addpath('../Libraries/cvx_a64');
% end
% 
% %% Add cvx to path
% if ~exist('cvx_begin', 'file')
%     cvx_setup;
% end

%% Add Peter Corke's toolbox to library path if required
if ~exist('rtbdemo', 'file')
    startup_rvc;
end