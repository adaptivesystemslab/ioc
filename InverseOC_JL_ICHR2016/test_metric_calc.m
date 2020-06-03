function test_metric_calc
% simulation metric calculation
clearvars
load('C:\Documents\MATLABResults\IOCProject\2016-03-17-15-15-47\sim.mat');

% pull out the weights and contribution ratios
len_weights = size(ccost_array{1}, 2);
repmat_fct = @(x) repmat(x, len_weights, 1);
sum_fct = @(x) sum(x);
max_fct = @(x) max(x);

c_doc_weight = cellfun(repmat_fct, ccost_array, 'UniformOutput', false);
weight_doc_full = vertcat(c_doc_weight{:});

c_doc_contrib = cellfun(repmat_fct, J_optsim_contrib_array, 'UniformOutput', false);
contrib_doc_full = vertcat(c_doc_contrib{:});

% pull out: rmse, norm, pivot, sign, condA, svdA, corrA 
rmse_ioc = [];
resnorm_ioc = [];
pivot_ioc = [];
condA_ioc = [];
corrVal_ioc = [];
svd_s_ioc = [];
weight_ioc = [];
contrib_ioc = [];
for ind_iter = 1:length(c_doc_weight)
    currOutput = output_inverse{ind_iter};
    currRmse = [currOutput(:).rmse]';
    currResnorm = [currOutput(:).resnorm]';
    
    [~, currMinRmseInd] = min(currRmse);
    [~, currMinResnormInd] = min(currResnorm);
    
    minRmseInd(ind_iter) = currMinRmseInd + (ind_iter-1)*len_weights;
    minResnormInd(ind_iter) = currMinResnormInd + (ind_iter-1)*len_weights;
    
    pivotSignSum = cellfun(sum_fct, {currOutput(:).gradDotSign}, 'UniformOutput', false);
    maxCorr = cellfun(max_fct, {currOutput(:).corrValArray}, 'UniformOutput', false);
    
    rmse_ioc = [rmse_ioc; [currOutput(:).rmse]'];
    resnorm_ioc = [resnorm_ioc; [currOutput(:).resnorm]'];
    condA_ioc = [condA_ioc; [currOutput(:).condA]'];
    pivot_ioc = [pivot_ioc; [pivotSignSum{:}]'];    
    
    corrVal_ioc = [corrVal_ioc; [maxCorr{:}]'];
    svd_s_ioc = [svd_s_ioc; [currOutput(:).svd_s]];
    
    weight_ioc = [weight_ioc; vertcat(currOutput(:).c_offset)];
    contrib_ioc = [contrib_ioc; vertcat(currOutput(:).J_out_contrib)];
end

c_ioc_weight_error = errorCalc(weight_doc_full, weight_ioc);
c_ioc_ratio_error = sum(errorCalc(contrib_doc_full, contrib_ioc), 2);

% determine the correlation between error and metric
h_main = figure; 
ax = subplot(231); 
generateAndPlot(ax, rmse_ioc, c_ioc_ratio_error, minRmseInd, minResnormInd, 'rmse');

ax = subplot(232); 
generateAndPlot(ax, resnorm_ioc, c_ioc_ratio_error, minRmseInd, minResnormInd, 'resnorm');

ax = subplot(233); 
generateAndPlot(ax, condA_ioc, c_ioc_ratio_error, minRmseInd, minResnormInd, 'cond');

ax = subplot(234); 
generateAndPlot(ax, pivot_ioc, c_ioc_ratio_error, minRmseInd, minResnormInd, 'pivot');

ax = subplot(235); 
generateAndPlot(ax, corrVal_ioc, c_ioc_ratio_error, minRmseInd, minResnormInd, 'max(corr)');

ax = subplot(236); 

ylabel('Contrib ratio err sum');
xlabel('Metric');

for ind_windowCount = 1:length(C_report_array)
    C_report_array{ind_windowCount}
%     gradOutput{i_la}
end
end

function result = errorCalc(base, target)
%     result = (base - target) ./ base;
    result = abs(base - target);
    
    result(isnan(result)) = 0;
    result(isinf(result)) = 0;
end

