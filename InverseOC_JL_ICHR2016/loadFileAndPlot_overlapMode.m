function [t_obs, avgWeightArray_ioc, avgRatioArray_ioc, cost_function_names_sorted, outputStruct] = loadFileAndPlot(matFilesToLoad, currInstName, indTotal, outputPathLocal, masterOverall, masterSave1)

individualplotflag = 0;
plotIOCWeightsOnly = 0;
plotUnnormedJflag = 1;
loadOnly=1;

    colorVec = [0,0.5,0; 0,0.8,0; 0.7,0.95,0; 1,0.8,0.1; 1,0.55,0; ...
            0.9,0,0.55; 0.4,0.1,0.7; 0.8,0.7,1; 0,0.5,0.9; 0,0.85,1; ...
            0.5,0.85,0.45; 0.7,0.6,0.2; 0.9,0.7,0.5; 0.8,0.4,0.3; 0.6,0.5,0.7; ...
            1,1,0; 1,0,1; 0,1,1; 0,0,1];
%     colorVec2 = [0,0.45,0.75; 0.85,0.3,0.1; 0.9,0.7,0.1; ...
%             0.5,0.2,0.55; 0.35,0.8,0.2; 0.4,0.75,0.95];
    colorVec2 = [0,0.5,0; 0,0.8,0; 0.7,0.95,0; 1,0.8,0.1; 1,0.55,0; ...
            0.9,0,0.55; 0.4,0.1,0.7; 0.8,0.7,1; 0,0.5,0.9; 0,0.85,1; ...
            0.5,0.85,0.45; 0.7,0.6,0.2; 0.9,0.7,0.5; 0.8,0.4,0.3; 0.6,0.5,0.7; ...
            1,1,0; 1,0,1; 0,1,1; 0,0,1];
x0=100;    y0=100;    width=900;    height=800;
    jointNames = {'q1', 'q2', 'q3', 'q4', 'q5', 'q6', 'q7'};
    jointsToPlot = [1:7]; % back_jFB and right side arm and leg joints

    outputPath = outputPathLocal;
    checkMkdir(outputPath);
    rmse_fct = @(a,b,c) sqrt(sum(sum((a-b).^2))/(size(a, 1)*size(a, 2)));

    thresholdMultiplier = 1;
    nowStr = datestr(now, 'yyyy_mm_dd_HH_MM_SS');

    plotCalc_assembleFile;
    plotCalc_metricCalc;
    
    [cost_function_names_sorted, cost_function_names_sorted_ind] = sort(cost_function_names);
    
    if individualplotflag
    for ind_indiv = 1:length(feature_win_save)
        t_obs = feature_win_save{ind_indiv}.t;
        q_obs = feature_win_save{ind_indiv}.q;
        dq_obs = feature_win_save{ind_indiv}.dq;
        ddq_obs = feature_win_save{ind_indiv}.ddq;
        
        t_doc = t_recon_plot_array{ind_indiv};
        q_doc = q_recon_plot_array{ind_indiv};
        dq_doc = dq_recon_plot_array{ind_indiv};
        ddq_doc = ddq_recon_plot_array{ind_indiv};
        
        weightSingle = output_inverse{ind_indiv}(minRmseIndArray(ind_indiv)).c_recovered;
        weightSingle = mean(weightSingle, 1);  weightSingle = weightSingle/sum(weightSingle);
        avgWeightArray_ioc = repmat(weightSingle, size(t_obs, 2), 1);
          
        ratioSingle = output_inverse{ind_indiv}(minRmseIndArray(ind_indiv)).J_out_contrib;
        ratioSingle = mean(ratioSingle, 1);  ratioSingle = ratioSingle/sum(ratioSingle);
        avgRatioArray_ioc = repmat(ratioSingle, size(t_obs, 2), 1);
        
        resnormSingle = output_inverse{ind_indiv}(minRmseIndArray(ind_indiv)).resnorm_lsqlin;
        resnromAll_lioc = repmat(resnormSingle, size(t_obs, 2), 1);
        
        rmsEntryToUseCurr = ones(size(t_obs));
        rmseMean = rmse_fct(q_obs, q_doc, []);
        rmseStd = 0;
        resnormMean = output_inverse{ind_indiv}(minRmseIndArray(ind_indiv)).resnorm_lsqlin;
        resnormStd = 0;
        
        % make individual plots
        outputStr = [pad(num2str(ind_indiv), 4, 'left', '0')];
        
        plot_series();
    end
    end
    
    % q_doc = q_blended
    % dq_doc = dq_blended
    % ddq_doc = ddq_blended
    
    % temp
    elapsedTime = elapsedTimeTotal;
    rmse_report.t = feature_full.t;
    rmse_report.maxminRMSE = [];
    rmse_report.meanRMSE = [];
    rmse_report.rangeRMSE = [];
    rmse_report.stdRMSE = [];
    t_doc = horzcat(t_recon_plot_array{:});
    q_doc = horzcat(q_recon_plot_array{:});
    dq_doc = horzcat(dq_recon_plot_array{:});
    ddq_doc = horzcat(ddq_recon_plot_array{:});
    
    t_obs = feature_full.t;
    q_obs = feature_full.q;
    dq_obs = feature_full.dq;
    ddq_obs = feature_full.ddq;
    
    dqtau_recon_plot_merge = horzcat(dqtau_recon_plot_array{:})';
    ddx_recon_plot_merge = horzcat(ddx_recon_plot_array{:})';
%     dqtau_plot_merge = horzcat(dqtau_plot_array{:})';
 

    
%     save(fullfile(outputPath, [currInstName '_postProcPackage']));
    
%     segmentInfo.timeStart = param.jump.takeoffFrame*0.005;
%     segmentInfo.timeEnd = param.jump.landFrame*0.005;
    
    % make individual plots
    avgWeightArray_ioc = avgWeightArray_belowThres;
    avgRatioArray_ioc = avgRatioArray_belowThres;
    resnromAll_lioc = resnromAll_lsqlin_const_minRMSE_array_belowThres;
    rmseMean = rmse_report.windowed_rmse_belowThreshold_mean;
    rmseStd = rmse_report.windowed_rmse_belowThreshold_std;
    rmsEntryToUseCurr = rmsEntryToUse;
    
    if(~individualplotflag) % need to declare these variables
        resnormMean = 0;
        resnormStd = 0;
    end
    
    outputStr = [pad(num2str(0000), 4, 'left', '0')];
    
    if ~(loadOnly)
        if(plotIOCWeightsOnly)
            plot_series_IOC_weights_only();
        else
            plot_series();
        end
    end
    
    outputStruct = struct;
%end
    
    % pull out the final window and normalize based on time
    timeStruct = createTaskRepetitionTimeStruct();
    ind = find(ismember([timeStruct.name], currInstName));
    dt = 0.01;
    
    % reconstruct the position data
    matDataPath = ['../../expressiveiocData/dataMat/2019_04_11_rightarm3/matEkfIk/' currInstName '_OnlyRightArm_' currInstName(8:9) '_mocap_mocap_X00_floating_ekfId1_ekfIk.mat'];
    modelPath = '../../../kalmanfilter/ik_framework/instance_expressiveioc/model/ioc_v4_rightarm_fixedbase.xml';
    modelInstance = rlModelInstance_expressiveioc_rightArm(0);
    modelInstance.loadModelFromModelSpecsNoSensor(modelPath, matDataPath);
    
    mdl = modelInstance.model;
%     mdl = rlCModel();
    f = mdl.getFrameByName('frame_rhand_end');
    pickId = round(timeStruct(ind).pickTarget/2);
    placementId = round(timeStruct(ind).placementTarget/2);
    indsToRun = pickId:placementId;
%     indsToRun = 1:length(feature_full.t);
    indX = 0;
    
    mdl.position = feature_full.q(:, indsToRun(1));
    mdl.forwardPosition();
    cart0 = f.t(1:3, 4);
    
    
%     if 0 
%         vis = rlVisualizer('vis',640,480);
%         mdl.forwardPosition();
%         vis.addModel(mdl);
%         vis.update();
%     end
    
    for ind_t = indsToRun
        indX = indX + 1;
        mdl.position = feature_full.q(:, ind_t);
        mdl.forwardPosition();
        cart = f.t(1:3, 4) - cart0;
        t_pos(indX) = feature_full.t(ind_t);
        x(indX,:) = cart;
        x_norm(indX) = norm(cart);
        
%         if 0
%             vis.update();
%         end
    end
    
    dx_norm = calcDeriv(x_norm, dt);
%     A = find(abs(dx_norm) < 0.05);
%     B = diff(A);
%     ii = [0, diff(B(:)')==0,0];
%     i1 = strfind(ii,[0 1]);
%     i2 = strfind(ii,[1 0]);
%     out = [B(i1)',i1(:),i2(:),i2(:)-i1(:)];
%     [~, longRange] = max(out(:,4));
%     midVal = floor((i2(longRange) - i1(longRange))/2);
%     midVal2 = midVal + sum(B(1:midVal)) - 2; % fix the diff offsets
    
    % midVal2 = timeStruct(ind).final.middle - timeStruct(ind).final.start;
    midVal2 = placementId - pickId;
    % midVal3 = timeStruct(ind).final.start + midVal2;
    
% %     if 0 
%     figure; subplot(211); plot(indsToRun, x_norm); hold on; 
%     plot(indsToRun(midVal2), x_norm(midVal2), '--gs',...
%         'LineWidth',2,...
%         'MarkerSize',10,...
%         'MarkerEdgeColor','b',...
%         'MarkerFaceColor',[0.5,0.5,0.5]);
%     title(currInstName);
%     subplot(212); plot(indsToRun, dx_norm); hold on;
%     plot(indsToRun(midVal2), dx_norm(midVal2), '--gs',...
%         'LineWidth',2,...
%         'MarkerSize',10,...
%         'MarkerEdgeColor','b',...
%         'MarkerFaceColor',[0.5,0.5,0.5]);
%     title(num2str(midVal3));
% %     end
    
    % now scale by position
    startingDist = 0;
    endingDist = x_norm(midVal2);
    x_norm_going = x_norm / endingDist;  startOffset = indsToRun(1);

    inc = 0:0.05:1;
    y = x_norm_going(1:midVal2)';
    x = linspace( 1, 0, length(y))';
    p = polyfit(x,y,1);
    pVal = polyval(p,x);
   
    indLa = [];
    weightLa = [];
    ratioLa = [];
    for i = 1:length(inc)
        errVal = abs(pVal - inc(i));
        [minErr, minInd] = min(errVal);
        
        indLa(i) = minInd;
%         [val, indLa(i)] = findClosestValue(inc(i), x_norm_going);
        ratioLa(:, i) = avgRatioArray_belowThres(startOffset + indLa(i), :);
        weightLa(:, i) = avgWeightArray_belowThres(startOffset + indLa(i), :);
    end
    
    outputStruct.currInstName = currInstName;
    outputStruct.indices_forward = indLa;
    outputStruct.ratio_forward = ratioLa;
    outputStruct.weights_forward = weightLa;

  %     Backward: placement to pick location (It is not needed any more)  
% %     inc = 1:-0.05:0;
% %     y = x_norm_going(midVal2:end)';
% %     x = linspace( 1, 0, length(y))';
% %     p = polyfit(x,y,1);
% %     pVal = polyval(p,x);
% %     
% %     indLa = [];
% %     weightLa = [];
% %     ratioLa = [];
% %     for i = 1:length(inc)
% %         errVal = abs(pVal - inc(i));
% %         [minErr, minInd] = min(errVal);
% %        
% %         indLa(i) = midVal2 + minInd;
% % %         [val, indLa(i)] = findClosestValue(inc(i), x_norm_going);
% %         ratioLa(:, i) = avgRatioArray_belowThres(startOffset + indLa(i), :);
% %         weightLa(:, i) = avgWeightArray_belowThres(startOffset + indLa(i), :);
% %     end
% %     
% %     outputStruct.indices_back = indLa;
% %     outputStruct.weights_back = weightLa;
% %     outputStruct.ratio_back = ratioLa;
    
%     diffRange = diff(lowVeloRange);
%     out = max(diff(find([1,diff(diffRange),1])));
    
    
%     if(plotUnnormedJflag)
%         % Plot "un-normalized" J weighting
%         if exist('h15', 'var') && h15 > 0
%             figure(h15);
%         else
%             h15 = figure;
%         end
% 
%         barObj = bar(feature_full.t(1:size(avgAbsArray,1)), avgAbsArray(:, cost_function_names_sorted_ind), 'stacked');
%         shading flat
%         for i = 1:numel(cost_function_names_sorted_ind)
%             barObj(i).FaceColor = colorVec(i,:);
%         end
%         yMaxHeight = max(sum(avgAbsArray,2));
%         ylim([0 yMaxHeight])
%         plotBoxes(gcf, segmentInfo.timeStart, segmentInfo.timeEnd, 'k', 0, yMaxHeight, -0.05);
%         ylabel('Recovered J absolute','Interpreter','Latex');
%         xlabel('Time [s]');
%         legend(cost_function_names_sorted,'Location','northwest');
%     %     xlim([0 8]);
%         grid on;
% 
%         set(gcf,'position',[x0,y0,width,height]);
%         saveas(h15, fullfile(outputPath, [currInstName '_J_abs_' outputStr '.fig']));
%         saveas(h15, fullfile(outputPath, [currInstName '_J_abs_' outputStr '.png']));
%     end
    
%     close all;
    
%     writeToOverallCSV(masterOverall, matFilesToLoad{1}, currInstName, elapsedTime, runSettings, rmse_report, cost_function_names, param);