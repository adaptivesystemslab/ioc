function IOCAnalysis()
    setPaths();
    
    perturbAmount = 1e-3;
    residualThreshold = 1e-3; % 0.1%

    
    nowstr = datestr(now, 'yyyymmdd_HHMMSS');
    basePath = 'D:\aslab_gitlab\expressive-ioc\Data\IOC\';
    masterPathCsv = ['D:\aslab_gitlab\expressive-ioc\Data\IOC\a_summary_' nowstr '.csv'];
    
    % iterate through all the folders in bathpath to produce figures and
    % table summary of the results
    
    % look for 'weights*1_' since that will be unique
    % then
    
    dirPath = dir(basePath);
    
    for i = 1:length(dirPath)
        currPath = fullfile(basePath, dirPath(i).name);
        
        if strcmpi(currPath(end), '.')
            continue; % it's . or ..
        end
        
        dirCurrSubPath = dir(currPath);
        
        for j = 1:length(dirCurrSubPath)
            currSubPath = fullfile(currPath, dirCurrSubPath(j).name);
            if strcmpi(currSubPath(end), '.')
                continue; % it's . or ..
            end
            
            if ~exist(currSubPath, 'dir')
                continue;
            end
            
%             if ~strcmpi(dirCurrSubPath(j).name, 'Subj1_3DOF_8CF')
%                 continue;
%             end
            
            matData = assembleData(currSubPath);
            
            suffix = [dirPath(i).name '_' matData.trialInfo.runName '_' matData.trialInfo.templateName];

%             try
                % plot results
                fprintf('%s\n', suffix);
                outputPathFig1 = fullfile(basePath, ['fig_results_individual_' suffix]);
                outputPathFig2 = fullfile(basePath, ['fig_results_cumulativeAllPass_' suffix]);
                outputPathFig3 = fullfile(basePath, ['fig_results_cumulativeRankPass_' suffix]);
                outputPathCsv = fullfile(basePath, ['csv_' suffix]);
                csv_populate(matData, masterPathCsv);
                plotting_individual(matData, outputPathFig1, outputPathCsv, masterPathCsv);
                plotting_cumulative(matData, outputPathFig2, outputPathFig3, outputPathCsv, masterPathCsv);
                
                % load data and perturb
%                 outputPathFig1 = fullfile(basePath, ['fig_perturb_single_' suffix]);
%                 outputPathFig2 = fullfile(basePath, ['fig_perturb_combined_' suffix]);
% %                 peturbWeights_single(matData, outputPathFig1, perturbAmount, residualThreshold);
%                 peturbWeights_combined(matData, outputPathFig2, perturbAmount, residualThreshold);
                
%             catch err
%                 err
%                 close all;
%             end
        end
    end
end

function peturbWeights_single(matData, outputPathFig1, perturbAmount, residualThreshold)
    progressVar = matData.progress;
    progressVarSecondary = matData.processSecondaryVar;
    featureLabels = matData.featureLabels;
    lenWeights = length(matData.featureLabels);
    t = matData.t;
    
    % remove only a single weight type at a time
    for i = 1:length(progressVar)
        weights(i, :) = progressVar(i).weights;
        
        H = progressVarSecondary(i).H;
        Hhat = H/norm(H,'fro');
        [weights_temp, residual_temp, x] = computeWeights(Hhat, lenWeights);
        
        % run baseline residual
        n=size(Hhat,2);
        H=Hhat'*Hhat;
        f=zeros(n,1);
        residual_baseline = 0.5*x'*H*x + f'*x;
        
        % perterb individual weights
        for j = 1:lenWeights
            x_perturb = x;
            x_perturb(j) = x_perturb(j) + perturbAmount;
            residual_peturb = 0.5*x_perturb'*H*x_perturb + f'*x_perturb;
            
%             residual_base(i) = residual_baseline;
%             residual_pert(i, :) = residual_peturb;
            residual_percent = abs(residual_peturb - residual_baseline)/residual_baseline;
            
            if residual_percent < residualThreshold
                % drop the entry
                residual_keep_single(i, j) = 0;
                x_perturb(j) = 0;
                x_perturb = x_perturb(1:lenWeights) / sum(x_perturb(1:lenWeights));
                weights_modified_single(i, :, j) = x_perturb;
            else
                residual_keep_single(i, j) = 1;
                weights_modified_single(i, :, j) = x(1:lenWeights);
            end
        end
    end
    
    % reassemble into progress var
    [weights_cum_orig, winCount_orig] = cumWeights(t, progressVar, featureLabels);
     
    for j = 1:lenWeights
        progressVarModified = progressVar;
        
        for i = 1:length(progressVarModified)
            progressVarModified(i).weights = weights_modified_single(i, :, j);
            progressVarModified(i).rankTraj = [];
            progressVarModified(i).rankPass = [];
            progressVarModified(i).error = [];
        end
        
         [weights_cum_modified{j}, winCount_cum_modified{j}] = cumWeights(t, progressVarModified, featureLabels);
    end
    
    h1 = figure('Position', [-1287 278 560*2 420*2]);
    ax(1) = subplot(lenWeights+1,2,1);
    area(t, weights);
    title('Individual original weights');
    
    ax(2) = subplot(lenWeights+1,2,2);
    area(t, weights_cum_orig);
    title('Cumulative original weights');
    
    for j = 1:lenWeights
        removeInds = find(residual_keep_single(:, j) == 0);
        
        axInd = (2*(j+1)-1);
        ax(axInd) = subplot(lenWeights+1,2,axInd);
        area_ax = area(t, weights_modified_single(:, :,j));
        hold on;
        plot_ax = plot(t(removeInds), -0.1*ones(size(removeInds)), 'LineStyle', 'none', 'Marker', '.', 'Color', area_ax(j).FaceColor);
        ylim([-0.2 1]);
        title(['Individual weights removing ' featureLabels{j} ', ' ...
            num2str(length(removeInds)) ' of ' num2str(size(weights_modified_single, 1)) ' removed']);
        
        axInd = (2*(j+1));
        ax(axInd) = subplot(lenWeights+1,2,axInd);
        area(t, weights_cum_modified{j});
        hold on;
        title(['Cumulative weights removing ' featureLabels{j} ', ' ...
            num2str(length(removeInds)) ' of ' num2str(size(weights_modified_single, 1)) ' removed']);
    end
    
    linkaxes(ax, 'x');
    xlim([0 t(end)]);
    
    saveas(h1, outputPathFig1, 'fig');
    saveas(h1, outputPathFig1, 'png');
    close(h1);
end 

function peturbWeights_combined(matData, outputPathFig2, perturbAmount, residualThreshold)
    progressVar = matData.progress;
    progressVarSecondary = matData.processSecondaryVar;
    featureLabels = matData.featureLabels;
    lenWeights = length(matData.featureLabels);
    t = matData.t;
    q = matData.q;
    dt = matData.dt;
    
    % now remove all the weights, starting from the smallest magnitude one
    for i = 1:length(progressVar)
        weights(i, :) = progressVar(i).weights;
        
        H = progressVarSecondary(i).H;
        Hhat = H/norm(H,'fro');
        [weights_temp, residual_temp, x] = computeWeights(Hhat, lenWeights);
        
        % run baseline residual
        n=size(Hhat,2);
        H=Hhat'*Hhat;
        f=zeros(n,1);
        residual_baseline = 0.5*x'*H*x + f'*x;
        
         % determine the smallest weight
        x_weights = x(1:lenWeights);
        [B, I] = sort(x_weights);
       
        for j_ind = 1:length(I)
            j = I(j_ind);
            x_perturb = x;
            x_perturb(j) = x_perturb(j) + perturbAmount;
            residual_peturb = 0.5*x_perturb'*H*x_perturb + f'*x_perturb;
            residual_percent = abs(residual_peturb - residual_baseline)/residual_baseline;
            
            if residual_percent < residualThreshold
                % drop the entry
                residual_keep_cumulative(i, j) = 0;
                x(j) = 0;
%                 x(1:lenWeights) = x(1:lenWeights) / sum(x(1:lenWeights));
            else
                residual_keep_cumulative(i, j) = 1;
            end
            
            residual_percent_cumulative(i, j) = residual_percent;
        end
        
        weights_modified_cumulative(i, :) = x(1:lenWeights) / sum(x(1:lenWeights));
    end
    
    % reassemble into progress var
    [weights_cum_orig, winCount_orig] = cumWeights(t, progressVar, featureLabels);
    
    progressVarModified = progressVar;
    for i = 1:length(progressVarModified)
        progressVarModified(i).weights = weights_modified_cumulative(i, :);
        progressVarModified(i).rankTraj = [];
        progressVarModified(i).rankPass = [];
        progressVarModified(i).error = [];
    end
    [weights_cum_cum, winCount_cum] = cumWeights(t, progressVarModified, featureLabels);
    
    totalCf = sum(residual_keep_cumulative, 1);
    totalKeep = sum(residual_keep_cumulative, 2);
    entriesZeros = find(totalKeep == 0);
    
    h2 = figure('Position', [-1287 278 560*2 420*2]);
    ax(1) = subplot(4,2,1);
    area_ax = area(t, weights);
    ylim([0 1]);
    title('Individual original weights');
    
    ax(2) = subplot(4,2,2);
    area(t, weights_cum_orig);
    ylim([0 1]);
    title('Cumulative original weights');
    
    ax(3) = subplot(4,2,3);
    area(t, weights_modified_cumulative); hold on;
    plot(t(entriesZeros), -0.1*ones(size(entriesZeros)), 'LineStyle', 'none', 'Marker', '.', 'Color', 'k');
    ylim([-0.2 1]);
    title(['Individual multiple-reduced weights, ' num2str(length(entriesZeros)) ' of ' num2str(numel(t)) ' all below thres']);
    
    ax(4) = subplot(4,2,4);
    area(t, weights_cum_cum);
    ylim([0 1]);
    title('Cumulative multiple-reduced weights');
    
    ax(5) = subplot(4,2,5);
%     hold on;
%     for i = 1:lenWeights
%         keptInds = find(residual_keep_cumulative(:, i) == 1);
%         plot_ax = plot(t(keptInds), i*ones(size(keptInds)), 'LineStyle', 'none', 'Marker', '.', 'Color', area_ax(i).FaceColor);
%     end
    area(t, residual_keep_cumulative);
    ylim([0 lenWeights+1]);
    title('Entries Kept');
    
    ax(6) = subplot(4,2,6);
    plot(t, q);
    
    ax(7) =  subplot(4,2,7); hold on
    for i = 1:lenWeights
        keptInds = find(residual_keep_cumulative(:, i) == 1);
        plot_ax = plot(t(keptInds), i*ones(size(keptInds)), 'LineStyle', 'none', 'Marker', '.', 'Color', area_ax(i).FaceColor);
    end
%     area(t, residual_keep_cumulative);
    ylim([0 lenWeights+1]);
    title('Entries Kept');
     
    ax(8) = subplot(4,2,8); hold on
%     c = categorical(featureLabels);
    for i = 1:length(featureLabels)
        barh(categorical(featureLabels(i)), totalCf(i)*dt, 'FaceColor', area_ax(i).FaceColor);
    end
%     line([numel(t) numel(t)], [0 length(featureLabels)]);

    linkaxes(ax, 'x');
    xlim([0 t(end)]);

    saveas(h2, outputPathFig2, 'fig');
    saveas(h2, outputPathFig2, 'png');
    close(h2);
end

function matData = assembleData(currSubPath)
    dirInnerPath = dir(fullfile(currSubPath, 'weights*.mat'));
    
    for i = 1:length(dirInnerPath)
        currFileName = dirInnerPath(i).name;
        datPath = fullfile(currSubPath, currFileName);
        load(datPath);
        currInds = outputVar_weights.frameInds(1):outputVar_weights.frameInds(end);
        matData.timeElapsed(i) = outputVar_weights.timeElapsed;
        matData.progress(currInds) = outputVar_weights.progress;
        
        currFileName = strrep(dirInnerPath(i).name, 'weights', 'data');
        datPath = fullfile(currSubPath, currFileName);
        load(datPath);
        matData.t(currInds) = outputVar_data.t;
        matData.q(currInds, :) = outputVar_data.q;
%         matData.dq(currInds, :) = outputVar_data.dq;
%         matData.tau(currInds, :) = outputVar_data.tau;
%         matData.features(currInds, :) = outputVar_data.features;
%         matData.dynamics(currInds, :) = outputVar_data.dynamics;

        currFileName = strrep(dirInnerPath(i).name, 'weights', 'supp');
        datPath = fullfile(currSubPath, currFileName);
        load(datPath);
        matData.processSecondaryVar(currInds) = outputVar_supp.processSecondaryVar;
    end
    
    matData.trialInfo = outputVar_weights.trialInfo;
    matData.timeElapsed = max(matData.timeElapsed);
    matData.minLenThres = outputVar_weights.minLenThres;
    matData.maxLenThres = outputVar_weights.maxLenThres;
    matData.minRankThres = outputVar_weights.minRankThres;
    
    matData.lenDof = outputVar_data.lenDof;
    matData.lenState = outputVar_data.lenState;
    matData.lenControl = outputVar_data.lenControl;
    matData.featureLabels = outputVar_data.featureLabels;
    matData.dt = outputVar_data.dt;
end

function csv_populate(matData, masterPathCsv)
    if ~exist(masterPathCsv, 'file')
        %             if doesn't exist write header
        header = 'Subject,Model,CF,Direction,DataFile';
        header = [header ',timeElapsedTotal,numFrames,timeElapsedPerLength'];
        header = [header '\n'];
    else
        header = '';
    end

    [fileID, errmsg] = fopen(masterPathCsv, 'a');
    fprintf(fileID, header);

    % write data to log file
    runNameSplit = strsplit(matData.trialInfo.runName, '_');
    subpath = strsplit(matData.trialInfo.subpath, '/');
    fprintf(fileID,'%s,%s,%s,%s,%s',runNameSplit{1}, runNameSplit{2}, ...
        runNameSplit{3}, matData.trialInfo.hWinAdvFlag, subpath{end});
    
%     framesProc = (matData.frameInds(end) - matData.frameInds(1) + 1);
    framesProc = length(matData.progress);
    timeElapsedPerFrame = matData.timeElapsed / framesProc;
    fprintf(fileID,',%f,%f,%f',matData.timeElapsed, framesProc, timeElapsedPerFrame);
    
    fprintf(fileID,'\n');
    fclose(fileID);
end

function h = plotting_cumulative(matData, outputPathFig_all, outputPathFig_rank, outputPathCsv, masterPathCsv)
    % load and process data
    progressVar = matData.progress;
%     dataIndsRan = matData.frameInds;
    rankThres = matData.minRankThres;
    featureLabels = matData.featureLabels;
    
%     if length(dataIndsRan) ~= length(progressVar)
        dataIndsRan = 1:length(progressVar);
%     end
    
    t = matData.t(dataIndsRan);
    q = matData.q(dataIndsRan, :);
    
    % the weight at each timestep is the sum of every window that overlaps with iterate
    [weights_all, winCount_all] = cumWeights(t, progressVar, featureLabels);
    [weights_rank, winCount_rank] = cumWeights_rankFilt(t, progressVar, featureLabels, rankThres);
    
    weightLabels = matData.featureLabels;
    
    % make fig
    h = figure('Position', [488 342 560*2 420*2]); 
    ax(1) = subplot(311);
    plot(t, q);
    ylabel('Joint Angle [rad]');
    title('Averaged weight windows (All pass)');
    
    ax(2) = subplot(312);
    area(t, weights_all);
    lgd = legend(weightLabels,'AutoUpdate','off');
    lgd.NumColumns = 1;
    ylabel('Weights [0:1]');
    
    ax(3) = subplot(313);
    area(t, winCount_all);
    ylabel('WinCountUsed');
    
    linkaxes(ax, 'x');
    xlabel('Time [s]');
    
    saveas(h, outputPathFig_all, 'fig');
    saveas(h, outputPathFig_all, 'png');
    close(h);
    
    
        % make fig
    h = figure('Position', [488 342 560*2 420*2]); 
    ax(1) = subplot(311);
    plot(t, q);
    ylabel('Joint Angle [rad]');
    title('Averaged weight windows (Rank threshold)');
    
    ax(2) = subplot(312);
    area(t, weights_rank);
    lgd = legend(weightLabels,'AutoUpdate','off');
    lgd.NumColumns = 1;
    ylabel('Weights [0:1]');
    
    ax(3) = subplot(313);
    area(t, winCount_rank);
    ylabel('WinCountUsed');
    
    linkaxes(ax, 'x');
    xlabel('Time [s]');
    
    saveas(h, outputPathFig_rank, 'fig');
    saveas(h, outputPathFig_rank, 'png');
    close(h);
end

function [weights_all, winCount_all] = cumWeights(t, progressVar, featureLabels)
    weights_all = zeros(length(t), length(featureLabels));
    winCount_all = zeros(length(t), 1);
    
    % the weight at each timestep is the sum of every window that overlaps with iterate
    for i = 1:length(t)
        weightAtI = [];
        for j = 1:length(progressVar)
            if progressVar(j).winInds(1) <= i && progressVar(j).winInds(end) >= i
                if sum(isnan(progressVar(j).weights)) == 0
                    weightAtI = [weightAtI; progressVar(j).weights];
                end
            end
        end
        
        weights_all(i, :) = mean(weightAtI, 1);
        winCount_all(i, :) = size(weightAtI, 1);
    end
end

function [weights_rank, winCount_rank] = cumWeights_rankFilt(t, progressVar, featureLabels, rankThres)
    weights_rank = zeros(length(t), length(featureLabels));
    winCount_rank = zeros(length(t), 1);
    
    for i = 1:length(t)
        weightAtI = [];
        for j = 1:length(progressVar)
            if progressVar(j).winInds(1) <= i && progressVar(j).winInds(end) >= i && progressVar(j).rankTraj(end) >= rankThres
                weightAtI = [weightAtI; progressVar(j).weights];
            end
        end
        
        if size(weightAtI, 1) > 0
            weights_rank(i, :) = mean(weightAtI, 1);
            winCount_rank(i, :) = size(weightAtI, 1);
        else
            weights_rank(i, :) = zeros(1, length(matData.featureLabels));
            winCount_rank(i, :) = 0;
        end
    end
end

function h = plotting_individual(matData, outputPathFig, outputPathCsv, masterPathCsv)
    % load and process data
    progressVar = matData.progress;
    minLenThres = matData.minLenThres;
    maxLenThres = matData.maxLenThres;
    minRankThres = matData.minRankThres;
    
%     dataIndsRan = matData.frameInds;
    
%     if length(dataIndsRan) ~= length(progressVar)
        dataIndsRan = 1:length(progressVar);
%     end
    
    t = matData.t(dataIndsRan);
    dt = matData.dt;
    q = matData.q(dataIndsRan, :);
    
    for i = 1:length(progressVar)
        weights(i, :) = progressVar(i).weights;
        rank(i, :) = progressVar(i).rankTraj(end);
        winLen(i, :) = progressVar(i).winInds(end) - progressVar(i).winInds(1) + 1;
    end
    
    weightLabels = matData.featureLabels;
    
    % make fig
    h = figure('Position', [488 342 560*2 420*2]); 
    ax(1) = subplot(411);
    plot(t, q);
    ylabel('Joint Angle [rad]');
    title('Individual weight windows');
    
    ax(2) = subplot(412);
    area(t, weights);
    lgd = legend(weightLabels,'AutoUpdate','off');
    lgd.NumColumns = 1;
    ylabel('Weights [0:1]');
    
    ax(3) = subplot(413);
    area(t, rank);
    hold on;
    plot([t(1) t(end)], [minRankThres minRankThres], 'k', 'LineWidth', 2);
    ylim([0 minRankThres*1.2]);
    ylabel('Final Rank');
    
    ax(4) = subplot(414);
    area(t, winLen*dt);
    hold on;
    plot([t(1) t(end)], [minLenThres minLenThres]*dt, 'k', 'LineWidth', 2);
    plot([t(1) t(end)], [maxLenThres maxLenThres]*dt, 'k', 'LineWidth', 2);
    ylim([minLenThres/1.2*dt maxLenThres*1.2*dt]);
    ylabel('Win Len [s]');
    
    linkaxes(ax, 'x');
    xlabel('Time [s]');
    
    saveas(h, outputPathFig, 'fig');
    saveas(h, outputPathFig, 'png');
    close(h);
end

function h = plotting(matData, outputPathFig, outputPathCsv, masterPathCsv)
    % load and process data
    t = matData.t;
    q = matData.q;
    weights = matData.weightTraj;
    rankRatio = matData.rankTraj;
    segArray = matData.segmentArray(1:end-1, :);
    weightLabels = matData.featureLabels;
    gamma = matData.runParam.gamma;
    
    if ~isempty(segArray)
        for i = 1:size(segArray, 1)
            segmentInfo(i).timeStart = t(segArray(i, 1));
            segmentInfo(i).timeEnd = t(segArray(i, 2));
            segmentInfo(i).timeWidth = segmentInfo(i).timeEnd - segmentInfo(i).timeStart;
            
            tTimeSeries(i) = segmentInfo(i).timeStart;
            weightTimeSeries(i, :) = weights(segArray(i, 2), :);
        end
        
        ind = size(segArray, 1);
        tTimeSeries(ind + 1) = segmentInfo(ind).timeEnd;
        weightTimeSeries(ind + 1, :) = weights(segArray(ind, 2), :);
        
        timeWidth = [segmentInfo.timeWidth];
    else
%         tTimeSeries(1) = 0;
%         tTimeSeries(2) = t(end);
%         weightTimeSeries(1, :) = zeros(size(weights, 2), 1);
%         weightTimeSeries(2, :) = zeros(size(weights, 2), 1);
    end
   
    % make fig
    h = figure('Position', [488 342 560*2 420*2]); 
    ax(1) = subplot(411);
    plot(t, q);
    if ~isempty(segArray)
        plotBoxes(segmentInfo, 'k');
    end
    ylabel('Joint Angle [rad]');
    
    ax(2) = subplot(412);
    if ~isempty(segArray)
        area(tTimeSeries, weightTimeSeries);
        lgd = legend(weightLabels,'AutoUpdate','off');
        lgd.NumColumns = 1;
        plotBoxes(segmentInfo, 'k', 0, -0.05, 1.05);
    end
    ylabel('Weights [0:1]');
    
    ax(3) = subplot(413);
    plot(t, weights);
    lgd = legend(weightLabels,'AutoUpdate','off');
    lgd.NumColumns = 1;
    if ~isempty(segArray)
        plotBoxes(segmentInfo, 'k', 0, -0.05, 1.05);
    end
    ylabel('Weights [0:1]');
    
    ax(4) = subplot(414);
    plot(t, rankRatio); 
    if ~isempty(segArray)
        ylim([-1 gamma]);
        plotBoxes(segmentInfo, 'k', 0, -10, gamma+10);
    end
    ylabel('Rank Ratio');
    
    linkaxes(ax, 'x');
    xlabel('Time [s]');
    
    saveas(h, outputPathFig, 'fig');
    saveas(h, outputPathFig, 'png');
    close(h);
end
