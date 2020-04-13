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