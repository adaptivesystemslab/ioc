function loadFileAndPlot(matFilesToLoad, indTotal, outputPathLocal, masterOverall, masterSave1)

    outputPath = outputPathLocal;
    checkMkdir(outputPath);
    rmse_fct = @(a,b,c) sqrt(sum(sum((a-b).^2))/(size(a, 1)*size(a, 2)));

    nowStr = datestr(now, 'yyyy_mm_dd_HH_MM_SS');
    nowStr = '';

    for ind = 1:length(matFilesToLoad)
        % load each individual file
        load(matFilesToLoad{ind});

    plotType = 'weight';
    plotCalc_metricCalc;
    
    % temp
    elapsedTime = 0;
    rmse_report.t = feature_full.t;
    rmse_report.maxminRMSE = max(minRmseVal);
    rmse_report.meanRMSE = mean(minRmseVal);
    rmse_report.rangeRMSE = range(minRmseVal);
    rmse_report.stdRMSE = std(minRmseVal);
    t_recon_plot_merge = horzcat(t_recon_plot_array{:})';
    q_recon_plot_merge = horzcat(q_recon_plot_array{:})';
    q_opt_plot_merge = feature_full.q;

%     writeToOverallCSV(masterOverall, currInstName, elapsedTime, runSettings, rmse_report, cost_function_names, param);
    
    plot_traj_weight;
    plot_resnorm_beforeThres;
    plot_resnorm_afterThres;
    plot_resnorm_bothThres;
    plot_paper_win;
    end
    
    saveas(h1, fullfile(outputPath, [currInstName '_fig1a_traj_' nowStr '.fig']));
    saveas(h1, fullfile(outputPath, [currInstName '_fig1a_traj_' nowStr '.png']));
    saveas(h11_1, fullfile(outputPath, [currInstName '_fig11_traj_' nowStr '_before.fig']));
    saveas(h11_1, fullfile(outputPath, [currInstName '_fig11_traj_' nowStr '_before.png']));
    saveas(h11_1, [masterSave1 '_before.fig']);
    saveas(h11_1, [masterSave1 '_before.png']);
    saveas(h11_2, fullfile(outputPath, [currInstName '_fig11_traj_' nowStr '_after.fig']));
    saveas(h11_2, fullfile(outputPath, [currInstName '_fig11_traj_' nowStr '_after.png']));
    saveas(h11_2, [masterSave1 '_after.fig']);
    saveas(h11_2, [masterSave1 '_after.png']);
    saveas(h11_3, fullfile(outputPath, [currInstName '_fig11_traj_' nowStr '_both.fig']));
    saveas(h11_3, fullfile(outputPath, [currInstName '_fig11_traj_' nowStr '_both.png']));
    saveas(h11_3, [masterSave1 '_after.fig']);
    saveas(h11_3, [masterSave1 '_after.png']);
    
    saveas(h14, fullfile(outputPath, [currInstName '_fig14_traj_' nowStr '.fig']));
    saveas(h14, fullfile(outputPath, [currInstName '_fig14_traj_' nowStr '.png']));
    saveas(h14, [masterSave1 '_winpaper.fig']);
    saveas(h14, [masterSave1 '_winpaper.png']);
    
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