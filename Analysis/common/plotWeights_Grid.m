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
            area_ax_pos = [left bottom width height];
            
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