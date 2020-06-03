% Only plot joint trajectories and recovered weights, no resnorm subplots
% and no q, dq, ddq plots.
% ALSO, saves figures in parent folder, not separate folder for each recording 

tmp = strsplit(outputPath,'\');
outputPathParent = strjoin(tmp(1:end-1),'\');

h1 = plot_fct_c_IOC_weights_only(feature_full, segmentInfo, ...
    t_obs, t_doc, q_obs, q_doc, avgWeightArray_ioc, resnromAll_lioc, ...
    rmsEntryToUseCurr, jointNames, jointsToPlot, colorVec, colorVec2, cost_function_names_sorted_ind, cost_function_names_sorted, ...
    rmseMean, rmseStd, resnormMean, resnormStd);
set(gcf,'position',[x0,y0,width,height]);
suptitle(['Jump Recording: ' strrep(name(1:9),'_','-')]);
saveas(h1, fullfile(outputPathParent, [currInstName '_c_' outputStr '.fig']));
saveas(h1, fullfile(outputPathParent, [currInstName '_c_' outputStr '.png']));

h2 = plot_fct_J_IOC_weights_only(feature_full, segmentInfo, ...
    t_obs, t_doc, q_obs, q_doc, avgRatioArray_ioc, resnromAll_lioc, ...
    rmsEntryToUseCurr, jointNames, jointsToPlot, colorVec, colorVec2, cost_function_names_sorted_ind, cost_function_names_sorted, ...
    rmseMean, rmseStd, resnormMean, resnormStd);
set(gcf,'position',[x0,y0,width,height]);
suptitle(['Jump Recording: ' strrep(name(1:9),'_','-')]);
saveas(h2, fullfile(outputPathParent, [currInstName '_J_' outputStr '.fig']));
saveas(h2, fullfile(outputPathParent, [currInstName '_J_' outputStr '.png']));

close all;