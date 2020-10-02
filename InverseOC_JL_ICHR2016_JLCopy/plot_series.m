keepInds = find(t_obs > 0);
t_obs = t_obs(keepInds);
q_obs = q_obs(:, keepInds);
dq_obs = dq_obs(:, keepInds);
ddq_obs = ddq_obs(:, keepInds);

avgWeightArray_ioc = avgWeightArray_ioc(keepInds, :);
avgWeightArray_belowThres = avgWeightArray_belowThres(keepInds, :);
resnromAll_lioc = resnromAll_lioc(keepInds);
avgRatioArray_ioc = avgRatioArray_ioc(keepInds, :);

h1 = plot_fct_c(feature_full, segmentInfo, ...
    t_obs, t_doc, q_obs, q_doc, avgWeightArray_ioc, resnromAll_lioc, ...
    rmsEntryToUseCurr, jointNames, jointsToPlot, colorVec, colorVec2, cost_function_names_sorted_ind, cost_function_names_sorted, ...
    rmseMean, rmseStd, resnormMean, resnormStd);
set(gcf,'position',[x0,y0,width,height]);
title(name);
saveas(h1, fullfile(outputPath, [currInstName '_c_' outputStr '.fig']));
saveas(h1, fullfile(outputPath, [currInstName '_c_' outputStr '.png']));

h2 = plot_fct_J(feature_full, segmentInfo, ...
    t_obs, t_doc, q_obs, q_doc, avgRatioArray_ioc, resnromAll_lioc, ...
    rmsEntryToUseCurr, jointNames, jointsToPlot, colorVec, colorVec2, cost_function_names_sorted_ind, cost_function_names_sorted, ...
    rmseMean, rmseStd, resnormMean, resnormStd);
set(gcf,'position',[x0,y0,width,height]);
title(name);
saveas(h2, fullfile(outputPath, [currInstName '_J_' outputStr '.fig']));
saveas(h2, fullfile(outputPath, [currInstName '_J_' outputStr '.png']));

h3 = plot_fct_q(feature_full, segmentInfo, ...
    t_obs, t_doc, q_obs, q_doc, dq_obs, dq_doc, ddq_obs, ddq_doc, ...
    rmsEntryToUseCurr, jointNames, jointsToPlot, colorVec, colorVec2, cost_function_names_sorted_ind, cost_function_names_sorted, ...
    rmseMean, rmseStd, resnormMean, resnormStd);
set(gcf,'position',[x0,y0,width,height]);
title(name);
saveas(h3, fullfile(outputPath, [currInstName '_q_' outputStr '.fig']));
saveas(h3, fullfile(outputPath, [currInstName '_q_' outputStr '.png']));

h4 = plot_fct_q_all(feature_full, segmentInfo, ...
    t_obs, t_doc, q_obs, q_doc, dq_obs, dq_doc, ddq_obs, ddq_doc, ...
    rmsEntryToUseCurr, jointNames, jointsToPlot, colorVec, colorVec2, cost_function_names_sorted_ind, cost_function_names_sorted, ...
    rmseMean, rmseStd, resnormMean, resnormStd);
set(gcf,'position',[x0,y0,width,height]);
% title(name);
saveas(h4, fullfile(outputPath, [currInstName '_q_all_' outputStr '.fig']));
saveas(h4, fullfile(outputPath, [currInstName '_q_all_' outputStr '.png']));

close all;