function loadFileAndPlot(matFilesToLoad, currInstName, indTotal, outputPathLocal, masterOverall, masterSave1)

    thresholdMultiplier = 1;

    outputPath = outputPathLocal;
    checkMkdir(outputPath);
    rmse_fct = @(a,b,c) sqrt(sum(sum((a-b).^2))/(size(a, 1)*size(a, 2)));

    nowStr = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
%     nowStr = '';
    
%     checkFile = fullfile(outputPath, [currInstName '_fig11_traj_' nowStr '_after.fig']);
%     if exist(checkFile, 'file')
%         return
%     end

    plotCalc_assembleFile;
    plotCalc_metricCalc;
    
    % temp
    elapsedTime = elapsedTimeTotal;
    rmse_report.t = feature_full.t;
    rmse_report.maxminRMSE = [];
    rmse_report.meanRMSE = [];
    rmse_report.rangeRMSE = [];
    rmse_report.stdRMSE = [];
    t_recon_plot_merge = horzcat(t_recon_plot_array{:})';
    q_recon_plot_merge = horzcat(q_recon_plot_array{:})';
    dq_recon_plot_merge = horzcat(dq_recon_plot_array{:})';
    ddq_recon_plot_merge = horzcat(ddq_recon_plot_array{:})';
    
    q_opt_plot_merge = feature_full.q;
    dq_opt_plot_merge = feature_full.dq;
    ddq_opt_plot_merge = feature_full.ddq;
    
    dqtau_recon_plot_merge = horzcat(dqtau_recon_plot_array{:})';
    ddx_recon_plot_merge = horzcat(ddx_recon_plot_array{:})';
%     dqtau_plot_merge = horzcat(dqtau_plot_array{:})';
 
    [cost_function_names_sorted, cost_function_names_sorted_ind] = sort(cost_function_names);
    
%     save(fullfile(outputPath, [currInstName '_postProcPackage']));
    
%     segmentInfo.timeStart = param.jump.takeoffFrame*0.005;
%     segmentInfo.timeEnd = param.jump.landFrame*0.005;
    
    colorVec = [0,0.5,0; 0,0.8,0; 0.7,0.95,0; 1,0.8,0.1; 1,0.55,0; ...
            0.9,0,0.55; 0.4,0.1,0.7; 0,0.5,1; 0,0.9,0.95; 0,1,0.9; 0.2,1,0.6];
    x0=100;    y0=100;    width=900;    height=800;
    
    plot_paper_win_3;
    set(gcf,'position',[x0,y0,width,height]);
    saveas(h14, fullfile(outputPath, [currInstName '_fig14q_c_' nowStr '.fig']));
    saveas(h14, fullfile(outputPath, [currInstName '_fig14q_c_' nowStr '.png']));
    
    plot_paper_3J; % modify me
    set(gcf,'position',[x0,y0,width,height]);
    saveas(h14, fullfile(outputPath, [currInstName '_fig14q_J_' nowStr '.fig']));
    saveas(h14, fullfile(outputPath, [currInstName '_fig14q_J_' nowStr '.png']));
    
%     plot_paper_sim;
%      saveas(h13, fullfile(outputPath, [currInstName '_fig1c_simc_' nowStr '.fig']));
%     saveas(h13, fullfile(outputPath, [currInstName '_fig1c_simc_' nowStr '.png']));
    
%     plot_paper_win_dqtau;
%         saveas(h14_dqtau, fullfile(outputPath, [currInstName '_fig14dqtau_traj_' nowStr '.fig']));
%     saveas(h14_dqtau, fullfile(outputPath, [currInstName '_fig14dqtau_traj_' nowStr '.png']));
%     
%     plot_paper_win_ddx;
%              saveas(h14_ddx, fullfile(outputPath, [currInstName '_fig14ddx_traj_' nowStr '.fig']));
%     saveas(h14_ddx, fullfile(outputPath, [currInstName '_fig14ddx_traj_' nowStr '.png']));
    
%     plot_correlation;
%     plot_correlation2;
    
%     plot_traj_weight;
%     set(gcf,'position',[x0,y0,width,height]);
%     saveas(h1, fullfile(outputPath, [currInstName '_fig1a_traj_' nowStr '.fig']));
%     saveas(h1, fullfile(outputPath, [currInstName '_fig1a_traj_' nowStr '.png']));
    
    
    % Plot "un-normalized" J weighting
    if exist('h15', 'var') && h15 > 0
        figure(h15);
    else
        h15 = figure;
    end

    barObj = bar(feature_full.t(1:size(avgAbsArray,1)), avgAbsArray(:, cost_function_names_sorted_ind), 'stacked');
    shading flat
    for i = 1:numel(cost_function_names_sorted_ind)
        barObj(i).FaceColor = colorVec(i,:);
    end
    yMaxHeight = max(sum(avgAbsArray,2));
    ylim([0 yMaxHeight])
    plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, yMaxHeight, -0.05);
    ylabel('Recovered $$\hat{J}$$ absolute','Interpreter','Latex');
    xlabel('Time [s]');
    legend(cost_function_names_sorted,'Location','northwest');
    xlim([0 8]);
    grid on;
    
    set(gcf,'position',[x0,y0,width,height]);
    saveas(h15, fullfile(outputPath, [currInstName '_fig14q_J_abs_' nowStr '.fig']));
    saveas(h15, fullfile(outputPath, [currInstName '_fig14q_J_abs_' nowStr '.png']));
    
    
    
    close all;
    
% % %     plot_trajWeight_individual;
    
   
%         plot_lambda;

    
%     plot_resnorm_beforeThres;
%     saveas(h11_1, fullfile(outputPath, [currInstName '_fig11_traj_' nowStr '_1before.fig']));
%     saveas(h11_1, fullfile(outputPath, [currInstName '_fig11_traj_' nowStr '_1before.png']));
%     
%     plot_resnorm_afterThres;
%     saveas(h11_2, fullfile(outputPath, [currInstName '_fig11_traj_' nowStr '_2after.fig']));
%     saveas(h11_2, fullfile(outputPath, [currInstName '_fig11_traj_' nowStr '_2after.png']));
    close all;
    
%     plot_paper_win;
  
   
%     plot_paper_win_dq;
%     plot_paper_win_ddq;
%     plot_ratio_standingvsmovement;
%     saveas(h4, fullfile(outputPath, [currInstName '_fig4a_plot_ratio_standingvsmovement_' nowStr '.fig']));
%     saveas(h4, fullfile(outputPath, [currInstName '_fig4a_plot_ratio_standingvsmovement_' nowStr '.png']));
    
%     plot_metrics2_small;
%     saveas(h2_2, fullfile(outputPath, [currInstName '_fig2b_line_' nowStr '.fig']));
%     saveas(h2_2, fullfile(outputPath, [currInstName '_fig2b_line_' nowStr '.png']));
%     saveas(h2_4, fullfile(outputPath, [currInstName '_fig2b_bar_' nowStr '.fig']));
%     saveas(h2_4, fullfile(outputPath, [currInstName '_fig2b_bar_' nowStr '.png']));
% 
%     plot_traj_weight_sliding;
%     saveas(h1_2, fullfile(outputPath, [currInstName '_fig1b_slidingtraj_' nowStr '.fig']));
%     saveas(h1_2, fullfile(outputPath, [currInstName '_fig1b_slidingtraj_' nowStr '.png']));
%     close all;
    
    if strcmpi(run_mode, 'sim')
        return
    end
    
    
    % by reps
% % %     if ~isfield(param, 'fatigue_level')
% % %         % no fatigue level defined
% % %         switch currFilestack.dataset
% % %             case 'squats_tuat_2011'
% % %                 switch currFilestack.subjectNumber
% % %                     case 1
% % %                         param.fatigue_level = [0 0 1 3 1];
% % %                         
% % %                     case 2
% % %                         param.fatigue_level = [2 1 4 5 1];
% % %                         
% % %                     case 3
% % %                         param.fatigue_level = [3 3 2 5 1];
% % %                         
% % %                     case 4
% % %                         param.fatigue_level = [2 2 3 5 1];
% % %                         
% % %                     case 5
% % %                         param.fatigue_level = [0 0 3 9 1];
% % %                         
% % %                     case 6
% % %                         param.fatigue_level = [2 3 7 12 1];
% % %                         
% % %                     case 7
% % %                         param.fatigue_level = [0 0 1 4 1];
% % %                         
% % %                     otherwise
% % %                         param.fatigue_level = [0 0 0 0 length(segmentInfo.timeStart)/5];
% % %                 end
% % %                 
% % %             otherwise
% % %                 param.fatigue_level = [0 0 0 0 length(segmentInfo.timeStart)/5];
% % %         end
% % %     elseif strcmpi(currFilestack.dataset, 'squats_tuat_2015')
% % %         param.fatigue_level = [0 0 0 0 length(segmentInfo.timeStart)/5];
% % % %     else
% % % %         param.fatigue_level = [0 0 0 0 length(segmentInfo.timeStart)];
% % %     end
% % % %     if sum(param.fatigue_level)*5 > length(segmentInfo.timeStart)
% % % %         % saved fatigue level is higher than the amount of segment defined
% % % %         param.fatigue_level = [0 0 0 0 length(segmentInfo.timeStart)/5];
% % % %     end
    
% % % if isfield(param, 'fatigue_level') && ~isempty(param.fatigue_level)
% % %     entriesToAdd = 1:((param.fatigue_level(1)-1) * 5) + 5;
% % %     if isempty(entriesToAdd)
% % %         entriesToAdd = 0;
% % %     end
% % %     plotGroups = {entriesToAdd};
% % %     for ind_rep = 2:length(param.fatigue_level)
% % %         entriesToAdd = plotGroups{ind_rep-1}(end) + [1:((param.fatigue_level(ind_rep)-1) * 5) + 5];
% % %         if isempty(entriesToAdd)
% % %             entriesToAdd = 0;
% % %         end
% % %         plotGroups{ind_rep} = entriesToAdd;
% % %     end
% % %     
% % %     try
% % %         [h18, rsq_mean_fatigue, rsq_std_fatigue] = plot_weightReps(figure, outputPath, [currInstName '_fig18_motionOnly_seg_'], feature_full, avgWeightArray_belowThres, resnromAll_lsqlin_const_minRMSE_array_belowThres, ...
% % %             cost_function_names, segmentInfo, plotGroups, 'motionOnly_seg');
% % %           title('Resnorm (motion), by fatigue');
% % %         
% % %         [h18l, rsq_mean_fatiguel, rsq_std_fatiguel] = plot_weightRepslambda(figure, outputPath, [currInstName '_fig18l_motionOnly_seg_'], feature_full, avgLambdaArray_belowThres, resnromAll_lsqlin_const_minRMSE_array_belowThres, ...
% % %             lambda_name, segmentInfo, plotGroups, 'motionOnly_seg');
% % %         
% % %         title('Resnorm (motion), by fatigue');
% % % %         h18_2 = plot_weightReps(figure, outputPath, [currInstName '_fig18_motionOnly_rest_'], feature_full, avgWeightArray_belowThres, resnromAll_lsqlin_const_minRMSE_array_belowThres, ...
% % % %             cost_function_names, segmentInfo, plotGroups, 'motionOnly_rest');
% % % %         title('Resnorm (rest), by fatigue');
% % %     catch err
% % %         fprintf(err.message);
% % %         h18 = figure;
% % %          h18l = figure;
% % %         rsq_mean_fatigue = zeros(size(cost_function_names, 2) + 1, 1);
% % %         rsq_std_fatigue = zeros(size(cost_function_names, 2) + 1, 1);
% % % %         h18_2 = figure;
% % %     end
% % %     saveas(h18, fullfile(outputPath, [currInstName '_fig18_traj_' nowStr '.fig']));
% % %     saveas(h18, fullfile(outputPath, [currInstName '_fig18_traj_' nowStr '.png']));
% % %     saveas(h18l, fullfile(outputPath, [currInstName '_fig18l_traj_' nowStr '.fig']));
% % %     saveas(h18l, fullfile(outputPath, [currInstName '_fig18l_traj_' nowStr '.png']));
% % %     
% % %     plotGroups = {};
% % %     segmentInfoJumps = 5;
% % %     for ind_rep = 1:segmentInfoJumps:length(segmentInfo.timeStart)
% % %         toAdd = (ind_rep):(ind_rep+segmentInfoJumps-1);
% % %         
% % %         if toAdd(end) > segmentInfoJumps:length(segmentInfo.timeStart)
% % %             toAdd = (ind_rep):length(segmentInfo.timeStart);
% % %         end
% % %         
% % %         plotGroups{end+1} = toAdd;
% % %     end
% % %     try
% % %         [h19, rsq_mean_set5, rsq_std_set5] = plot_weightReps(figure, outputPath, [currInstName '_fig19_motionOnly_seg_'], feature_full, avgWeightArray_belowThres, resnromAll_lsqlin_const_minRMSE_array_belowThres, ...
% % %             cost_function_names, segmentInfo, plotGroups, 'motionOnly_seg');
% % %         title('Resnorm (motion), by sets of 5');
% % %         
% % %         [h19l, rsq_mean_set5l, rsq_std_set5l] = plot_weightRepslambda(figure, outputPath, [currInstName '_fig19l_motionOnly_seg_'], feature_full, avgLambdaArray_belowThres, resnromAll_lsqlin_const_minRMSE_array_belowThres, ...
% % %             lambda_name, segmentInfo, plotGroups, 'motionOnly_seg');
% % %         title('Resnorm (motion), by sets of 5');        
% % % %         h19_2 = plot_weightReps(figure, outputPath, [currInstName '_fig19_motionOnly_rest_'], feature_full, avgWeightArray_belowThres, resnromAll_lsqlin_const_minRMSE_array_belowThres, ...
% % % %             cost_function_names, segmentInfo, plotGroups, 'motionOnly_rest');
% % % %         title('Resnorm (rest), by sets of 5');
% % %     catch err
% % %            fprintf(err.message);
% % %         h19 = figure;
% % %         h19l = figure;
% % %         rsq_mean_set5 = zeros(size(cost_function_names, 2) + 1, 1);
% % %         rsq_std_set5 = zeros(size(cost_function_names, 2) + 1, 1);
% % % %         h19_2 = figure;
% % %     end
% % %         saveas(h19, fullfile(outputPath, [currInstName '_fig19_traj_' nowStr '.fig']));
% % %     saveas(h19, fullfile(outputPath, [currInstName '_fig19_traj_' nowStr '.png']));
% % %             saveas(h19l, fullfile(outputPath, [currInstName '_fig19l_traj_' nowStr '.fig']));
% % %     saveas(h19l, fullfile(outputPath, [currInstName '_fig19l_traj_' nowStr '.png']));
% % %     close all;
% % % else
% % %     rsq_mean_set5 = zeros(size(cost_function_names, 2) + 10, 1);
% % %     rsq_std_set5 = zeros(size(cost_function_names, 2) + 10, 1);
% % %     rsq_mean_fatigue = zeros(size(cost_function_names, 2) + 10, 1);
% % %     rsq_std_fatigue = zeros(size(cost_function_names, 2) + 10, 1);
% % %     
% % %     rsq_mean_set5l = zeros(size(cost_function_names, 2) + 100, 1);
% % %     rsq_std_set5l = zeros(size(cost_function_names, 2) + 100, 1);
% % %     rsq_mean_fatiguel = zeros(size(cost_function_names, 2) + 100, 1);
% % %     rsq_std_fatiguel = zeros(size(cost_function_names, 2) + 100, 1);
% % % end

    rsq_mean_set5 = zeros(size(cost_function_names, 2) + 10, 1);
    rsq_std_set5 = zeros(size(cost_function_names, 2) + 10, 1);
    rsq_mean_fatigue = zeros(size(cost_function_names, 2) + 10, 1);
    rsq_std_fatigue = zeros(size(cost_function_names, 2) + 10, 1);
    
    rsq_mean_set5l = zeros(size(cost_function_names, 2) + 100, 1);
    rsq_std_set5l = zeros(size(cost_function_names, 2) + 100, 1);
    rsq_mean_fatiguel = zeros(size(cost_function_names, 2) + 100, 1);
    rsq_std_fatiguel = zeros(size(cost_function_names, 2) + 100, 1);
    
    

    rmse_report.rsq_mean_set5 = rsq_mean_set5;
    rmse_report.rsq_std_set5 = rsq_std_set5;
    rmse_report.rsq_mean_fatigue = rsq_mean_fatigue;
    rmse_report.rsq_std_fatigue = rsq_std_fatigue;
    
    rmse_report.rsq_namel = lambda_name;
    rmse_report.rsq_mean_set5l = rsq_mean_set5l;
    rmse_report.rsq_std_set5l = rsq_std_set5l;
    rmse_report.rsq_mean_fatiguel = rsq_mean_fatiguel;
    rmse_report.rsq_std_fatiguel = rsq_std_fatiguel;
    
    writeToOverallCSV(masterOverall, matFilesToLoad{1}, currInstName, elapsedTime, runSettings, rmse_report, cost_function_names, param);
    
%     h20 = plot_resnormrmse(output_inverse, minRmseIndArray, indToUse_window, currLoadPackage.minRmseVal);
%         saveas(h20, fullfile(outputPath, [currInstName '_fig20_traj_' nowStr '.fig']));
%     saveas(h20, fullfile(outputPath, [currInstName '_fig20_traj_' nowStr '.png']));
    
%     plot_metrics;
%     saveas(h2_1, fullfile(outputPath, [currInstName '_fig2a_metrics_' nowStr '.fig']));
%     saveas(h2_1, fullfile(outputPath, [currInstName '_fig2a_metrics_' nowStr '.png']));
%    
   
%     plot_metrics3;
%     saveas(h2_3, fullfile(outputPath, [currInstName '_fig2c_metrics_' nowStr '.fig']));
%     saveas(h2_3, fullfile(outputPath, [currInstName '_fig2c_metrics_' nowStr '.png']));
% 
%     
%     plot_metrics4;
%     saveas(h2_4, fullfile(outputPath, [currInstName '_fig2d_metrics_' nowStr '.fig']));
%     saveas(h2_4, fullfile(outputPath, [currInstName '_fig2d_metrics_' nowStr '.png']));
  
    
%     plot_individratios;
%     saveas(h3_1, fullfile(outputPath, [currInstName '_fig3a_individratios_' nowStr '.fig']));
%     saveas(h3_1, fullfile(outputPath, [currInstName '_fig3a_individratios_' nowStr '.png']));
%     close(h3_1);
%     
%     plot_individratios2;
%     saveas(h3_2, fullfile(outputPath, [currInstName '_fig3b_individratios_' nowStr '.fig']));
%     saveas(h3_2, fullfile(outputPath, [currInstName '_fig3b_individratios_' nowStr '.png']));
%     close(h3_2);
%     
%     plot_com_correlation;
%     saveas(h5_1, fullfile(outputPath, [currInstName '_fig5a_correlation_' nowStr '.fig']));
%     saveas(h5_1, fullfile(outputPath, [currInstName '_fig5a_correlation_' nowStr '.png']));
%     close(h5_1);
%     
%     plot_com_correlation2;
%     saveas(h5_2, fullfile(outputPath, [currInstName '_fig5b_correlation_' nowStr '.fig']));
%     saveas(h5_2, fullfile(outputPath, [currInstName '_fig5b_correlation_' nowStr '.png']));
%     close(h5_2);
%     
% % % % %     plot_ratio_standingvsmovement;
% % % % %     
% % % % % %     plot_recoveredweights;
% % % % % %     saveas(h6, fullfile(outputPath, [currInstName '_fig6a_traj_' datestr(now, 'yyyy_mm_dd_HH_MM_SS') '.fig']));
% % % % % %     saveas(h6, fullfile(outputPath, [currInstName '_fig6a_traj_' datestr(now, 'yyyy_mm_dd_HH_MM_SS') '.png']));
% % % % % %     close(h6);
% % % % %     
% % % % %     plot_resnorm12;
% % % % %     plot_minNegConst;
% % % % %     plot_resnormThreshold;
% % % % %     plot_resnorm_beforeThres;
% % % % %     plot_resnorm_afterThres;
% % % % %     plot_resnorm_bothThres;
% % %     plot_paper_sim;
% % % % %     plot_paper_win;
    
   
  

%     saveas(h7, fullfile(outputPath, [currInstName '_fig7a_traj_' nowStr '.fig']));
%     saveas(h7, fullfile(outputPath, [currInstName '_fig7a_traj_' nowStr '.png']));
%     saveas(h9, fullfile(outputPath, [currInstName '_fig9_traj_' nowStr '.fig']));
%     saveas(h9, fullfile(outputPath, [currInstName '_fig9_traj_' nowStr '.png']));
%     saveas(h10, fullfile(outputPath, [currInstName '_fig10_traj_' nowStr '.fig']));
%     saveas(h10, fullfile(outputPath, [currInstName '_fig10_traj_' nowStr '.png']));

%     saveas(h11_1, [masterSave1 '_before.fig']);
%     saveas(h11_1, [masterSave1 '_before.png']);

%     saveas(h11_2, [masterSave1 '_after.fig']);
%     saveas(h11_2, [masterSave1 '_after.png']);
%     saveas(h11_3, fullfile(outputPath, [currInstName '_fig11_traj_' nowStr '_both.fig']));
%     saveas(h11_3, fullfile(outputPath, [currInstName '_fig11_traj_' nowStr '_both.png']));
%     saveas(h11_3, [masterSave1 '_after.fig']);
%     saveas(h11_3, [masterSave1 '_after.png']);
%     saveas(h13, fullfile(outputPath, [currInstName '_fig13_traj_' nowStr '.fig']));
%     saveas(h13, fullfile(outputPath, [currInstName '_fig13_traj_' nowStr '.png']));
%     saveas(h13, [masterSave1 '_simpaper.fig']);
%     saveas(h13, [masterSave1 '_simpaper.png']);



%     saveas(h14_dq, fullfile(outputPath, [currInstName '_fig14dq_traj_' nowStr '.fig']));
%     saveas(h14_dq, fullfile(outputPath, [currInstName '_fig14dq_traj_' nowStr '.png']));
%     saveas(h14_ddq, fullfile(outputPath, [currInstName '_fig14ddq_traj_' nowStr '.fig']));
%     saveas(h14_ddq, fullfile(outputPath, [currInstName '_fig14ddq_traj_' nowStr '.png']));
%     saveas(h14, [masterSave1 '_winpaper.fig']);
%     saveas(h14, [masterSave1 '_winpaper.png']);
  

    
% % % % % %     saveas(h18_2, fullfile(outputPath, [currInstName '_fig18_2_traj_' nowStr '.fig']));
% % % % % %     saveas(h18_2, fullfile(outputPath, [currInstName '_fig18_2_traj_' nowStr '.png']));
% % % % % %     saveas(h19_2, fullfile(outputPath, [currInstName '_fig19_2_traj_' nowStr '.fig']));
% % % % % %     saveas(h19_2, fullfile(outputPath, [currInstName '_fig19_2_traj_' nowStr '.png']));
    

%     close(h7);
%     close(h9);    close(h10); close(h11_1); close(h11_2);  close(h11_3);     close(h14);     close(h4);    close(h2_4);  close(h2_3); close(h2_2);
%     close(h2_1);
    close all;
    
%     plot_stuff;
%     for ind_h8 = 1:length(h8)
%         saveas(h8(ind_h8), fullfile(outputPath, [currInstName '_fig8a_' cost_function_names{ind_h8} '_' nowStr '.fig']));
%         saveas(h8(ind_h8), fullfile(outputPath, [currInstName '_fig8a_' cost_function_names{ind_h8} '_' nowStr '.png']));
%         close(h8(ind_h8));      
%     end

%     plot_finalContribRatioWeight;
%     saveas(h12, fullfile(outputPath, [currInstName '_fig12_traj_' nowStr '_comp.fig']));
%     saveas(h12, fullfile(outputPath, [currInstName '_fig12_traj_' nowStr '_comp.png']));
%     saveas(h12, [masterSave1 '_comp.fig']);
%     saveas(h12, [masterSave1 '_comp.png']);
%     close(h12);