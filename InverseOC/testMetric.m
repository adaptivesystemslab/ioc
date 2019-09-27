%% set up environment
addpath(genpath('.'));

%% Load data
load('C:\Users\jf2lin\Downloads\wetransfer-a4cb56\result01\Subj1_3DOF_4CF_struct_weights.mat');
load('C:\Users\jf2lin\Downloads\wetransfer-a4cb56\result01\Subj1_3DOF_4CF_struct_data.mat');
load('C:\Users\jf2lin\Downloads\wetransfer-a4cb56\result01\Subj1_3DOF_4CF_struct_supp.mat');

outputVar_supp_4cf_data = outputVar_data;
outputVar_supp_4cf_weight = outputVar_weights;
outputVar_supp_4cf_supp = outputVar_supp;

load('C:\Users\jf2lin\Downloads\wetransfer-a4cb56\result01\Subj1_3DOF_8CF_struct_weights.mat');
load('C:\Users\jf2lin\Downloads\wetransfer-a4cb56\result01\Subj1_3DOF_8CF_struct_data.mat');
load('C:\Users\jf2lin\Downloads\wetransfer-a4cb56\result01\Subj1_3DOF_8CF_struct_supp.mat');

outputVar_supp_8cf_data = outputVar_data;
outputVar_supp_8cf_weight = outputVar_weights;
outputVar_supp_8cf_supp = outputVar_supp;

%% Look at the H matrix
for i = 1:5000
    H_4cf = genH(outputVar_supp_4cf_supp.processSecondaryVar(i).H);
    H_8cf = genH(outputVar_supp_8cf_supp.processSecondaryVar(i).H);

    [weights_temp4, residual_4cf(i, :), x_4cf(i, :), H_test_4cf] = computeWeights(H_4cf, 4);
    [weights_temp8, residual_8cf(i, :), x_8cf(i, :), H_test_8cf] = computeWeights(H_8cf, 8);
    
    weights_4cf(i, :) = weights_temp4';
    weights_8cf(i, :) = weights_temp8';
    
    range_4cf(i, :) = max(max(H_4cf)) - min(min(H_4cf));
    range_8cf(i, :) = max(max(H_8cf)) - min(min(H_8cf));
end

t = outputVar_supp_4cf_data.t(1:5000);

figure(1);
ax_weight(1) = subplot(4, 2, 1);
area(t, weights_4cf); ylabel('weights'); 
ax_residual(1) = subplot(4, 2, 3);
area(t, residual_4cf); ylabel('residual'); 
ax_x(1) = subplot(4, 2, 5); 
area(t, x_4cf); ylabel('all_x');
ax_range(1) = subplot(4, 2, 7);
plot(t, range_4cf); ylabel('range')

ax_weight(2) = subplot(4, 2, 2);
area(t, weights_8cf);  
ax_residual(2) = subplot(4, 2, 4);
area(t, residual_8cf);  
ax_x(2) = subplot(4, 2, 6); 
area(t, x_8cf); 
ax_range(2) = subplot(4, 2, 8);
plot(t, range_8cf);

linkaxes([ax_weight ax_residual ax_x ax_range], 'x');
% xlim([0 0.5]);

% H_4cf = H_4cf / min(min(abs(H_4cf)));
% H_8cf = H_8cf / min(min(abs(H_8cf)));

% compare range to weight



%% misc H
function Hhat = genH(H)
% H = outputVar_supp_4cf_supp.processSecondaryVar(1).H;
    Hhat = H/norm(H,'fro');
end