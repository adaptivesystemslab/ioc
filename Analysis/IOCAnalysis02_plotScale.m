function IOCAnalysis()
    setPaths();
    faceColours = brewermap(8, 'paired');
    
%     nowstr = datestr(now, 'yyyymmddHHMMSS');
    nowstr = '20200413_FatigueFull_3CF';
%     basePath = 'D:\aslab_gitlab\expressive-ioc\Data\IOC\';
    basePath = 'D:\results\fatigue_ioc01_weightsIndividual\20200413_FatigueFull_3CF\';
    outputPath = ['D:\results\fatigue_ioc02_weightsAssembled\' nowstr];
    masterPathCsv = [outputPath '\a_summary.csv'];
    checkMkdir(outputPath);
    
    currSubPathDir = dir(basePath);
    
        for j = 1:length(currSubPathDir)
            currSubPath = fullfile(basePath, currSubPathDir(j).name);
            if strcmpi(currSubPath(end), '.')
                continue; % it's . or ..
            end
            
            if ~exist(currSubPath, 'dir')
                continue;
            end
            
            matData = assembleData(currSubPath); 
%             suffix = [matData.trialInfo.runName '_' matData.trialInfo.templateName]
            suffix = [matData.trialInfo.runName];
            matData.t = (1:length(matData.t)) * 0.01;
            
            outputPathFig1 = fullfile(outputPath, ['fig_weiInd_' suffix]);
            outputPathFig2 = fullfile(outputPath, ['fig_weiCum_' suffix]);
            outputPathFig3 = fullfile(outputPath, ['fig_results_cumulativeRankPass_' suffix]);
            
            outputPathMat0 = fullfile(outputPath, ['mat_dataInd_' suffix]);
            outputPathMat1 = fullfile(outputPath, ['mat_weiCum_' suffix]);
            outputPathMat2 = fullfile(outputPath, ['mat_weiInd_' suffix]);
            
            if exist([outputPathMat1 '.mat'], 'file')
                fprintf('%s detected, skipping\n', outputPathMat1);
                continue;
            end            
            
            fprintf('Processing %s\n', outputPathMat1);
            
%             try
                % plot results
%                 fprintf('%s\n', suffix);
              
                save(outputPathMat0, 'matData');
                outputPathCsv = fullfile(outputPath, ['csv_' suffix]);
%                 csv_populate(matData, masterPathCsv);
                plotting_individual(matData, outputPathFig1, outputPathCsv, masterPathCsv, outputPathMat2);
                matSave = plotting_cumulative(matData, outputPathFig2, outputPathFig3, outputPathCsv, masterPathCsv, faceColours, outputPathMat1);
                
%                 % load data and perturb
%                 for ind_perturbAmount = 1:length(perturbAmount)
%                     for ind_residual = 1:length(residualThreshold)
%                         suffix_suffix = [num2str(ind_perturbAmount) '_' num2str(ind_residual)];
%                         suffix2 = [suffix '_' suffix_suffix];
%                         outputPathFig1 = fullfile(outputPath, ['fig_perturb_single_' suffix2]);
%                         outputPathFig2 = fullfile(outputPath, ['fig_perturb_combined_' suffix2]);
%                         suffix2 = [suffix '_' num2str(ind_perturbAmount) '_' num2str(ind_residual)]
%                         outputPathFig1 = fullfile(outputPath, ['fig_perturb_single_' suffix2]);
%                         outputPathFig2 = fullfile(outputPath, ['fig_perturb_combined_' suffix2]);
%                         matFilePath = fullfile(outputPath, ['mat_perturb_combined_' suffix2]);
%                         %                        peturbWeights_single(matData, outputPathFig1, currPerturb, currResidual);
%                         
%                         [t{ind_perturbAmount, ind_residual}, weights_cum_cum{ind_perturbAmount, ind_residual}, weights_blocks{ind_perturbAmount, ind_residual}, ...
%                             weight_labels{ind_perturbAmount, ind_residual}, residual_keep_cumulative{ind_perturbAmount, ind_residual}] = ...
%                             peturbWeights_combined(matData, outputPathFig2, perturbAmount(ind_perturbAmount), residualThreshold(ind_residual), faceColours, suffix_suffix, ...
%                             matSave, matFilePath);
%                     end
%                 end
                
%                 figFileName = fullfile(outputPath, ['a_fig_perturb_combinedGridWeight_' suffix]);
%                 plotWeights_Grid(figFileName, perturbAmount, residualThreshold, maxAcross, maxDown, faceColours,  ...
%                     t, weights_cum_cum, residual_keep_cumulative, [], weight_labels);
%                 
%                 figFileName = fullfile(outputPath, ['a_fig_perturb_combinedGridRemoval_' suffix]);
%                 plotWeights_Grid(figFileName, perturbAmount, residualThreshold, maxAcross, maxDown, faceColours,  ...
%                     t, [], residual_keep_cumulative, weights_blocks, weight_labels);
                       
%             catch err
%                 err
%                 close all;
%             end
        end
end

function plotWeights_Grid(figFileName, perturbAmount, residualThreshold, maxAcross, maxDown, faceColours,  ...
    t, weights_cum_cum, residual_keep_cumulative, weights_blocks, weight_labels)

    numRecsAcross = length(perturbAmount);
    numRecsDown = length(residualThreshold);
    winSepAcross = ceil(numRecsAcross/maxAcross);
    winSepDown = ceil(numRecsDown/maxDown);
    numFigs = (winSepAcross)*(winSepDown);
    
    % cue up the figures
    for ind_figNum = 1:numFigs
        h_grid(ind_figNum) = figure('Position', [22.6000 52.2000 1.8608e+03 924]);
    end
    
    % then assign the grid
    figMat = zeros(numRecsAcross, numRecsDown);
    currFigInd = 0;
    for ind_x = 1:winSepAcross
        for ind_y = 1:winSepDown
            currFigInd = currFigInd + 1;
            matInd_x = (ind_x-1)*maxAcross + [1:maxAcross];
            matInd_y = (ind_y-1)*maxDown + [1:maxDown];

            matInd_x = matInd_x(matInd_x <= numRecsAcross);
            matInd_y = matInd_y(matInd_y <= numRecsDown);
            figMat(matInd_x, matInd_y) = currFigInd;
        end
    end
    
    ind_overallInds = 0;
%     masterout = [];
    ax = [];
    for ind_perturbAmount = 1:numRecsAcross
        for ind_residual = 1:numRecsDown
            ind_overallInds = ind_overallInds + 1;
            ind_x = mod(ind_perturbAmount, maxAcross);
            if ind_x == 0
                ind_x = maxAcross;
            end
            ind_y = mod(ind_residual, maxDown);
            if ind_y == 0
                ind_y = maxDown;
            end
            currFig = figMat(ind_perturbAmount, ind_residual);
            currH = h_grid(currFig);
            figure(currH);
            
            left = (ind_x-1)*1/maxAcross;
            bottom = (maxDown-ind_y)*1/maxDown;
            width = (1/maxAcross) - 0.01;
            height = (1/maxDown) - 0.03;
            area_ax_pos = [left bottom width height];
%             masterout = [masterout; currFig ind_x ind_y area_ax_pos]
            ax(currFig, ind_overallInds) = subplot('Position', area_ax_pos);
            
            % insert plot
            if ~isempty(weights_cum_cum) 
                area_ax = area(t{ind_perturbAmount, ind_residual}, weights_cum_cum{ind_perturbAmount, ind_residual});
                ylim([0 1]);
                
                for ind_face = 1:length(area_ax)
                    area_ax(ind_face).FaceColor = faceColours(ind_face, :);
                end
            else
                for ind1 = 1:size(weights_blocks{ind_perturbAmount, ind_residual}, 3)
                    ack = weights_blocks{ind_perturbAmount, ind_residual}(:, :, ind1);
                    v = [];
                    f = [];
                    for ind2 = 1:size(ack, 1)
                        hold on;
                        currF = size(f, 1)*4;
                        if ack(ind2, 1) == 0
                            continue
                        end
                        x1 = t{ind_perturbAmount, ind_residual}(ack(ind2, 1));
                        x2 = t{ind_perturbAmount, ind_residual}(ack(ind2, 2));
                        y1 = ind1-1;
                        y2 = ind1;
                        v = [v; x1 y1; x1 y2; x2 y2; x2 y1];
                        f = [f; [1:4]+currF];
                    end
                    
                    patch('Faces', f, 'Vertices', v, 'FaceColor', faceColours(ind1, :));
                end
            end
            
            xlim([0 t{1}(end)]);
            
            totalKeep = sum(residual_keep_cumulative{ind_perturbAmount, ind_residual}, 2);
            entriesZeros = find(totalKeep == 0);
            title(['P ' num2str(perturbAmount(ind_perturbAmount)) ' | T ' num2str(residualThreshold(ind_residual)) ...
                ' | ABT ' num2str(length(entriesZeros)) '/' num2str(numel(t{1}))]);
            
            if ind_x == 1 && ind_y == 1
                legend(weight_labels{end},'AutoUpdate','off', 'Location','northwest');
            end
        end
    end
                
    for ind_figNum = 1:numFigs
%         linkaxes(ax(ind_figNum, :));
%         legend(weight_labels{end},'AutoUpdate','off');
        
        outputPathFig = [figFileName '_' num2str(ind_figNum)];
        saveas(h_grid(ind_figNum), outputPathFig, 'fig');
        saveas(h_grid(ind_figNum), outputPathFig, 'png');
        close(h_grid(ind_figNum));
    end
    
    % then make a master version of the plot
    h = figure('Position', [22.6000 52.2000 1.8608e+03 924]);
    ind_overallInds = 0;
     for ind_perturbAmount = 1:numRecsAcross
        for ind_residual = 1:numRecsDown
            ind_overallInds = ind_overallInds + 1;
            left = (ind_perturbAmount-1)*1/numRecsAcross;
            bottom = (numRecsDown-ind_residual)*1/numRecsDown;
            width = (1/numRecsAcross) - 0.005;
            height = (1/numRecsDown) - 0.005;
            area_ax_pos = [left bottom width height]
            
            ax(currFig, ind_overallInds) = subplot('Position', area_ax_pos);
            
             % insert plot
            if ~isempty(weights_cum_cum) 
                area_ax = area(t{ind_perturbAmount, ind_residual}, weights_cum_cum{ind_perturbAmount, ind_residual});
                ylim([0 1]);
                
                for ind_face = 1:length(area_ax)
                    area_ax(ind_face).FaceColor = faceColours(ind_face, :);
                end
            else
                for ind1 = 1:size(weights_blocks{ind_perturbAmount, ind_residual}, 3)
                    ack = weights_blocks{ind_perturbAmount, ind_residual}(:, :, ind1);
                    v = [];
                    f = [];
                    for ind2 = 1:size(ack, 1)
                        hold on;
                        currF = size(f, 1)*4;
                        if ack(ind2, 1) == 0
                            continue
                        end
                        x1 = t{ind_perturbAmount, ind_residual}(ack(ind2, 1));
                        x2 = t{ind_perturbAmount, ind_residual}(ack(ind2, 2));
                        y1 = ind1-1;
                        y2 = ind1;
                        v = [v; x1 y1; x1 y2; x2 y2; x2 y1];
                        f = [f; [1:4]+currF];
                    end
                    
                    patch('Faces', f, 'Vertices', v, 'FaceColor', faceColours(ind1, :));
                end
            end
            
            xlim([0 t{1}(end)]);
            
            totalKeep = sum(residual_keep_cumulative{ind_perturbAmount, ind_residual}, 2);
            entriesZeros = find(totalKeep == 0);
            title(['P ' num2str(perturbAmount(ind_perturbAmount)) ' | T ' num2str(residualThreshold(ind_residual)) ...
                ' | ABT ' num2str(length(entriesZeros)) '/' num2str(numel(t{1}))]);
        end
     end
     
     outputPathFig = [figFileName '_0'];
     saveas(h, outputPathFig, 'fig');
     saveas(h, outputPathFig, 'png');
     close(h);
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
    [weights_cum_orig, weights_cum_orig_var, winCount_orig] = cumWeights(t, progressVar, featureLabels);
     
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

function [t, weights_cum_cum, weights_blocks, featureLabels, residual_keep_cumulative] = peturbWeights_combined(matData, outputPathFig2, perturbAmount, residualThreshold, faceColours, suffix_suffix, matSave_unperturb, outputPathMat1)
    progressVar = matData.progress;
    progressVarSecondary = matData.processSecondaryVar;
    featureLabels = matData.featureLabels;
    lenWeights = length(matData.featureLabels);
    t = matData.t;
    q = matData.q;
    dt = matData.dt;
    weightLabels = matData.featureLabels;
    
    % now remove all the weights, starting from the smallest magnitude one
    for i = 1:length(progressVar)
        if isempty(progressVar(i).weights)
            weights(i, :) = zeros(size(progressVar(1).weights));
            residual_keep_cumulative(i, lenWeights) = 0;
            weights_modified_cumulative(i, :) = zeros(size(progressVar(1).weights));
            continue;
        end
        
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
    
    weights_blocks(1, 2, size(residual_keep_cumulative, 2)) = 0;
    for i = 1:size(residual_keep_cumulative, 2)
        finds = find(residual_keep_cumulative(:, i) == 1);
        
        if isempty(finds)
            continue
        end
        
        currVal = finds(1);
        weightBlockInds = 1;
        weights_blocks(weightBlockInds, 1, i) = currVal;
         
        for j = 2:length(finds)
            if finds(j) == currVal+1
                currVal = finds(j);
            else
                weights_blocks(weightBlockInds, 2, i) = currVal;
                
                weightBlockInds = weightBlockInds + 1;
                currVal = finds(j);
                weights_blocks(weightBlockInds, 1, i) = currVal;
            end
        end
        
        weights_blocks(weightBlockInds, 2, i) = finds(end);
     end
    
    % reassemble into progress var
%     [weights_cum_orig, weights_cum_var_orig, winCount_orig] = cumWeights(t, progressVar, featureLabels);
    weights_cum_orig = matSave_unperturb.weights;
    
    progressVarModified = progressVar;
    for i = 1:length(progressVarModified)
        progressVarModified(i).weights = weights_modified_cumulative(i, :);
        progressVarModified(i).rankTraj = [];
        progressVarModified(i).rankPass = [];
        progressVarModified(i).error = [];
    end
    [weights_cum_cum, weights_cum_cum_var, winCount_cum] = cumWeights(t, progressVarModified, featureLabels);
    
    matSave.t = t;
    matSave.q = q;
    matSave.weights = weights_cum_cum;
    matSave.var = weights_cum_cum_var;
    matSave.winCount = winCount_cum;
    save(outputPathMat1, 'matSave');
    
    
    totalCf = sum(residual_keep_cumulative, 1);
    totalKeep = sum(residual_keep_cumulative, 2);
    entriesZeros = find(totalKeep == 0);
    
    h2 = figure('Position', [-1287 278 560*2 420*2]);
    ax(1) = subplot(4,2,1);
    area_ax = area(t, weights);
    ylim([0 1]);
    title('Individual original weights');
    for i = 1:length(area_ax)
        area_ax(i).FaceColor = faceColours(i, :);
    end
    
    ax(2) = subplot(4,2,2);
    area_ax = area(t, weights_cum_orig);
    ylim([0 1]);
    title('Cumulative original weights');
    for i = 1:length(area_ax)
        area_ax(i).FaceColor = faceColours(i, :);
    end
    
    ax(3) = subplot(4,2,3);
    area_ax = area(t, weights_modified_cumulative); hold on;
    plot(t(entriesZeros), -0.1*ones(size(entriesZeros)), 'LineStyle', 'none', 'Marker', '.', 'Color', 'k');
    ylim([-0.2 1]);
    title(['Individual multiple-reduced weights, ' num2str(length(entriesZeros)) ' of ' num2str(numel(t)) ' all below thres']);
    for i = 1:length(area_ax)
        area_ax(i).FaceColor = faceColours(i, :);
    end
    
    ax(4) = subplot(4,2,4);
    area_ax = area(t, weights_cum_cum);
    ylim([0 1]);
    title('Cumulative multiple-reduced weights');
    for i = 1:length(area_ax)
        area_ax(i).FaceColor = faceColours(i, :);
    end
    
    ax(5) = subplot(4,2,6);
    plot(t, q);
    
    ax(6) = subplot(4,2,5); hold on
%     c = categorical(featureLabels);
    for i = 1:length(featureLabels)
        barh(categorical(featureLabels(i)), totalCf(i)*dt, 'FaceColor', area_ax(i).FaceColor);
    end

    ax(7) = subplot(4, 2, [7 8]);
    for i = 1:size(residual_keep_cumulative, 2)  
        ack = weights_blocks(:, :, i);
        v = [];
        f = [];
        for j = 1:size(ack, 1)
            hold on;
            currF = size(f, 1)*4;
            if ack(j, 1) == 0
                continue
            end
            x1 = t(ack(j, 1));
            x2 = t(ack(j, 2));
            y1 = i-1;
            y2 = i;
            v = [v; x1 y1; x1 y2; x2 y2; x2 y1];
            f = [f; [1:4]+currF];
        end
        
        patch('Faces', f, 'Vertices', v, 'FaceColor', area_ax(i).FaceColor);
    end
    title(['Entries Kept, perturb by ' num2str(perturbAmount) ', remove threshold on ' num2str(residualThreshold)]);

    linkaxes(ax, 'x');
    xlim([0 t(end)]);

    saveas(h2, outputPathFig2, 'fig');
    saveas(h2, outputPathFig2, 'png');
    close(h2);

    outputPathFig_all_variance = [outputPathFig2 '_var'];
    plot_weightVar(weights_cum_cum, weights_cum_cum_var, t,q,faceColours, outputPathFig_all_variance, weightLabels);
end

function plot_weightVar(weights_all, weights_all_var, t,q,faceColours, outputPathFig_all_variance, weightLabels)
 h = figure('Position', [488 342 560*2 420*2]);
    for i = 1:size(weights_all, 2)+1
        switch size(weights_all, 2)
            case 3
                ax(i) = subplot(2,2,i);
            case 4
                ax(i) = subplot(2,3,i);
            case 8
                ax(i) = subplot(3,3,i);
        end
        
        switch i
            case 1
                plot(t, q);
                title('q');
            otherwise
                ind_weight = i-1;
%                 boundedline(t,weights_all(:,ind_weight), weights_all_var(:,ind_weight), 'cmap', faceColours(ind_weight, :));
                boundedline(t,weights_all(:,ind_weight), weights_all_var(:,ind_weight), '-b');
                title(weightLabels{ind_weight});
                ylim([-0.25 1.25]);
        end
    end
    
    linkaxes(ax, 'x');
    xlabel('Time [s]');
    
    saveas(h, outputPathFig_all_variance, 'fig');
    saveas(h, outputPathFig_all_variance, 'png');
    close(h);
end

function matData = assembleData(currSubPath)
    iterPath = fullfile(currSubPath, 'weights*.mat');
    dirInnerPath = dir(iterPath);
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
        matData.dq(currInds, :) = outputVar_data.dq;
        matData.tau(currInds, :) = outputVar_data.tau;
        matData.features(currInds, :) = outputVar_data.features;
        matData.dynamics(currInds, :) = outputVar_data.dynamics;

        currFileName = strrep(dirInnerPath(i).name, 'weights', 'supp');
        datPath = fullfile(currSubPath, currFileName);
        load(datPath);
        matData.processSecondaryVar(currInds) = outputVar_supp.processSecondaryVar;
    end
    
    iterPath2 = strrep(iterPath, 'END', 'START');
    currSubPath = strrep(currSubPath, 'END', 'START');
    dirInnerPath = dir(iterPath2);
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
        matData.dq(currInds, :) = outputVar_data.dq;
        matData.tau(currInds, :) = outputVar_data.tau;
        matData.features(currInds, :) = outputVar_data.features;
        matData.dynamics(currInds, :) = outputVar_data.dynamics;

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

function [matSave] = plotting_cumulative(matData, outputPathFig_all, outputPathFig_rank, outputPathCsv, masterPathCsv, faceColours, outputPathMat1)
    outputPathFig_all_overall = [outputPathFig_all '_overall'];
    outputPathFig_all_variance = [outputPathFig_all '_variance'];
    
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
    [weights_all, weights_all_var, winCount_all] = cumWeights(t, progressVar, featureLabels);
%     [weights_rank, weights_rank_var, winCount_rank] = cumWeights_rankFilt(t, progressVar, featureLabels, rankThres);

    weightLabels = matData.featureLabels;

%     matSave = matData;
    matSave.t = t;
    matSave.q = q;
    matSave.weights = weights_all;
    matSave.var = weights_all_var;
    matSave.winCount = winCount_all;
    save(outputPathMat1, 'matSave');
    
    startInds = [1 100];
    endInds = max(t) - [100 0];
    
    % make fig
    h = figure('Position', [488 342 560*2 420*2]); 
    ax(1) = subplot(321);
    plot(t, q);
    ylabel('Joint Angle [rad]');
    title('Averaged weight windows (All pass)');
    xlim(startInds);
    
    ax(2) = subplot(322);
    plot(t, q);
    ylabel('Joint Angle [rad]');
    title('Averaged weight windows (All pass)');
    xlim(endInds);
    
    ax(3) = subplot(323);
    area(t, weights_all);
    lgd = legend(weightLabels,'AutoUpdate','off');
    lgd.NumColumns = 1;
    ylabel('Weights [0:1]');
    xlim(startInds);
    
    ax(4) = subplot(324);
    area(t, weights_all);
    lgd = legend(weightLabels,'AutoUpdate','off');
    lgd.NumColumns = 1;
    ylabel('Weights [0:1]');
    xlim(endInds);
    
    ax(5) = subplot(325);
    area(t, winCount_all);
    ylabel('WinCountUsed');
    xlim(startInds);
    
    ax(6) = subplot(326);
    area(t, winCount_all);
    ylabel('WinCountUsed');
    xlim(endInds);
%     linkaxes(ax, 'x');
%     xlabel('Time [s]');
    
    saveas(h, outputPathFig_all_overall, 'fig');
    saveas(h, outputPathFig_all_overall, 'png');
    close(h);
    
    % make fig
    plot_weightVar(weights_all, weights_all_var, t,q,faceColours, outputPathFig_all_variance, weightLabels);
     
    %         % make fig
    %     h = figure('Position', [488 342 560*2 420*2]);
    %     ax(1) = subplot(311);
    %     plot(t, q);
    %     ylabel('Joint Angle [rad]');
    %     title('Averaged weight windows (Rank threshold)');
    %
    %     ax(2) = subplot(312);
    %     area(t, weights_rank);
    %     lgd = legend(weightLabels,'AutoUpdate','off');
    %     lgd.NumColumns = 1;
    %     ylabel('Weights [0:1]');
    %
    %     ax(3) = subplot(313);
    %     area(t, winCount_rank);
    %     ylabel('WinCountUsed');
    %
%     linkaxes(ax, 'x');
%     xlabel('Time [s]');
%     
%     saveas(h, outputPathFig_rank, 'fig');
%     saveas(h, outputPathFig_rank, 'png');
%     close(h);
end

function [weights_mean_all, weights_var_all, winCount_all] = cumWeights(t, progressVar, featureLabels)
    weights_mean_all = zeros(length(t), length(featureLabels));
    weights_var_all = zeros(length(t), length(featureLabels));
    winCount_all = zeros(length(t), 1);
    
    % the weight at each timestep is the sum of every window that overlaps with iterate
    for i = 1:length(t)
        if mod(i, 1000) ==  0
            fprintf('[%u/%u] Processing weights\n', i, length(t));
        end
        
        weightAtI = [];
        for j = 1:length(progressVar)
            if isempty(progressVar(j).winInds)
                continue
            end
            
            if progressVar(j).winInds(1) <= i && progressVar(j).winInds(end) >= i
                if sum(isnan(progressVar(j).weights)) == 0
                    weightAtI = [weightAtI; progressVar(j).weights];
                end
            end
        end
        
        if ~isempty(weightAtI)
            weights_mean_all(i, :) = mean(weightAtI, 1);
            weights_var_all(i, :) = std(weightAtI, 0, 1);
            winCount_all(i, :) = size(weightAtI, 1);
        else
            weights_mean_all(i, :) = zeros(size(progressVar(1).weights));
            weights_var_all(i, :) = zeros(size(progressVar(1).weights));
            winCount_all(i, :) = 0;
        end
    end
end

function [weights_mean_rank, weights_var_rank, winCount_rank] = cumWeights_rankFilt(t, progressVar, featureLabels, rankThres)
    weights_mean_rank = zeros(length(t), length(featureLabels));
    weights_var_rank = zeros(length(t), length(featureLabels));
    winCount_rank = zeros(length(t), 1);
    
    for i = 1:length(t)
        weightAtI = [];
        for j = 1:length(progressVar)
            if isempty(progressVar(j).winInds)
                continue
            end
            
            if progressVar(j).winInds(1) <= i && progressVar(j).winInds(end) >= i && progressVar(j).rankTraj(end) >= rankThres
                weightAtI = [weightAtI; progressVar(j).weights];
            end
        end
        
        if size(weightAtI, 1) > 0
            weights_mean_rank(i, :) = mean(weightAtI, 1);
            weights_var_rank(i, :) = std(weightAtI, 1);
            winCount_rank(i, :) = size(weightAtI, 1);
        else
            weights_mean_rank(i, :) = zeros(1, length(featureLabels));
            winCount_rank(i, :) = 0;
        end
    end
end

function h = plotting_individual(matData, outputPathFig, outputPathCsv, masterPathCsv, outputPathMat2)
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
        if ~isempty(progressVar(i).weights)
            weights(i, :) = progressVar(i).weights;
            rank(i, :) = progressVar(i).rankTraj(end);
            winLen(i, :) = progressVar(i).winInds(end) - progressVar(i).winInds(1) + 1;
        else
            weights(i, :) = zeros(size(matData.featureLabels));
            rank(i, :) = 0;
            winLen(i, :) = 0;
        end
    end
    
    weightLabels = matData.featureLabels;
    
    % make fig
    h = figure('Position', [-1.4014e+03 153.8000 560*2 420*2]); 
%     ax(1) = subplot(311);
%     plot(t, q);
%     ylabel('Joint Angle [rad]');
%     title('Individual weight windows');
%     
%     ax(1) = subplot(211);
    area(t, weights);
    lgd = legend(weightLabels,'AutoUpdate','off');
    lgd.NumColumns = 1;
    ylabel('Weights [0:1]');
    ylim([0 1]);
    
%     ax(3) = subplot(413);
%     area(t, rank);
%     hold on;
%     plot([t(1) t(end)], [minRankThres minRankThres], 'k', 'LineWidth', 2);
%     ylim([0 minRankThres*1.2]);
%     ylabel('Final Rank');
    
%     ax(4) = subplot(414);
%     area(t, winLen*dt);
%     hold on;
%     plot([t(1) t(end)], [minLenThres minLenThres]*dt, 'k', 'LineWidth', 2);
%     plot([t(1) t(end)], [maxLenThres maxLenThres]*dt, 'k', 'LineWidth', 2);
%     ylim([minLenThres/1.2*dt maxLenThres*1.2*dt]);
%     ylabel('Win Len [s]');
    
    
    data = cell2mat({progressVar.winInds}');
    minimum = data(:, 1);
    maximum = data(:, 2);
%     h = figure;
    
%     ax(2) = subplot(212);
%     aH = axes;

%     bH = bar(t, maximum);
%     hold on;
%     bH1 = bar(t, minimum);
%     for ii = bH1
%         ii.FaceColor = [1 1 1];
%         ii.LineStyle = 'None';
%         ii.EdgeColor = [1 1 1];
%     end
%     for ii = bH
%         ii.LineStyle = 'None';
%         ii.EdgeColor = [1 1 1];
%     end
% %     axes('Position',aH.Position,'XTick',[],'YTick',[],'Color','None');
% 
%    linkaxes(ax, 'x');
%    xlim([0 t(end)]);
%     xlabel('Time [s]');
%     
    saveas(h, outputPathFig, 'fig');
    saveas(h, outputPathFig, 'png');
    close(h);
    
    matSave.t = t;
    matSave.q = q;
    matSave.weights = weights;
    save(outputPathMat2, 'matSave');
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
