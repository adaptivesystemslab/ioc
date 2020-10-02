% plot data
[~, uniqueInd] = unique(t_rmse);
[~, closeInd1] = findClosestValue(t_rmse(1), feature_full.t);
[~, closeInd2] = findClosestValue(t_rmse(end), feature_full.t);

h2_4 = figure;
ax = [];

q_toExamine = feature_full.q(1, closeInd1:closeInd2);
dq_toExamine = feature_full.dq(1, closeInd1:closeInd2);

q_lowerRange = min(q_toExamine);
q_upperRange = max(q_toExamine);
q_dRange = (q_upperRange - q_lowerRange)/10;

dq_lowerRange = -100;
dq_upperRange = 100;
dq_dRange = (dq_upperRange - dq_lowerRange)/2;

q_range = q_lowerRange:q_dRange:q_upperRange;
dq_range = dq_lowerRange:dq_dRange:dq_upperRange;

titles = {};
ratioMean = [];
ratioSD = [];
for ind_combos2 = 1:length(dq_range)-1
    for ind_combos1 = 1:length(q_range)-1
        q_lower = q_range(ind_combos1);
        q_upper = q_range(ind_combos1+1);
        dq_lower = dq_range(ind_combos2);
        dq_upper = dq_range(ind_combos2+1);
        titles = [titles ['q:' num2str(q_lower) ':' num2str(q_upper) ',dq:' num2str(dq_lower) ':' num2str(dq_upper)]];
        
        ind_q1 = (q_toExamine < q_upper);
        ind_q2 = (q_toExamine > q_lower);
        ind_dq1 = (dq_toExamine < dq_upper);
        ind_dq2 = (dq_toExamine > dq_lower);
        ind_int = find(ind_q1 + ind_q2 + ind_dq1 + ind_dq2 == 4);
        
        mean_val = mean(ratioRMSE_plot(ind_int, :));
        sd_val = std(ratioRMSE_plot(ind_int, :));
        
        if sum(isnan(mean_val))
            mean_val = zeros(size(mean_val));
            sd_val = zeros(size(sd_val));
        end
        
        ratioMean(end+1, :) = mean_val;
        ratioSD(end+1, :) = sd_val;
    end
end

% ind_pq = find(q_toExamine > 0);
% ind_pdq = find(dq_toExamine > 0);
% ind_nq = find(q_toExamine < 0);
% ind_ndq = find(dq_toExamine < 0);
% 
% titles = {'+q +dq', '+q -dq', '-q +dq', '-q -dq'};
% ind_interset{1} = intersect(ind_pq, ind_pdq); % positive q, positive dq
% ind_interset{2} = intersect(ind_pq, ind_ndq); % positive q, negative dq
% ind_interset{3} = intersect(ind_nq, ind_pdq); % negative q, positive dq
% ind_interset{4} = intersect(ind_nq, ind_ndq); % negative q, negative dq
% 
% for ind_combos = 1:length(titles)
%     ratioMean(ind_combos, :) = mean(ratioRMSE_plot(ind_interset{ind_combos}, :));
%     ratioSD(ind_combos, :) = std(ratioRMSE_plot(ind_interset{ind_combos}, :));
% end

h = barwitherr(ratioSD', ratioMean');% Plot with errorbars
legend(titles);
% xlim([0.5 length(cost_function_names)+0.5]);
set(gca,'XTickLabel',cost_function_names);
title('Contrib ratio in +-qdq');