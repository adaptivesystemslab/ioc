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
            weights(i, :) = zeros(1, lenWeights);
            residual_keep_cumulative(i, lenWeights) = 0;
            weights_modified_cumulative(i, :) = zeros(1, lenWeights);
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