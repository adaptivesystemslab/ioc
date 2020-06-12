function IOCAnalysis()
    setPaths();
    nownowstr = datestr(now, 'yyyymmddHHMMSS');
    sourceSuffix = '20200413_FatigueFull_3CF';
    targetSuffix = '20200413_FatigueFull_3CF_3';
    searchString = 'mat_*_3DOF_3CF*.mat';
      
    basePath = ['D:\results\fatigue_ioc03_weightsPattern\' sourceSuffix '\mat\'];
    outputPath = ['D:\results\fatigue_ioc04_weightsCluster\' targetSuffix '\'];
    
    outputPathStd = [outputPath 'std\'];
    outputPathReg = [outputPath 'reg\'];
    
    outStdCsv = [outputPath 'analysis_standard_' nownowstr '.csv'];
    outRegCsv = [outputPath 'analysis_regression_' nownowstr '.csv'];
    checkMkdir(outputPath);
    
    mean_mag_threshold = 1e-2;
    close all;
    
    currBasePathDir = dir([basePath searchString]);
    for j = 1:length(currBasePathDir)
        currFileName = currBasePathDir(j).name;
        currFullPath = fullfile(basePath, currFileName);
        
        if strcmpi(currFileName(end), '.')
            continue; % it's . or ..
        end

        % load each entry and prepare to plot
        load(currFullPath);
        
        subjectNum = str2num(trialInfo.runName(8:9));
%         combinedStats = [stats_q stats_dq stats_tau stats_dtau stats_weights stats_dweights];
        combinedStats = [stats_weights stats_dweights];
        cumStats{subjectNum} = combinedStats;
    end
    
    nSubject = length(cumStats);
    allDofs = {cumStats{2}.name};
    allFeaturesSingle = cumStats{2}(1).segStats_SingleWindow.Properties.VariableNames(5:end);
    
%     % pull correlation data
%     for ind_subject = 1:length(cumStats)
%         currSubj = cumStats{ind_subject};
%         
%         if isempty(currSubj)
%             continue
%         end
%         
%         for ind_dof = 1:length(allDofs)
% %             currStats = cumStats{ind_subject}(ind_dof);
%             currFeatureTable = cumStats{ind_subject}(ind_dof).segStats_SingleWindow;
%             currIndivTable = cumStats{ind_subject}(ind_dof).regression_individual_singleWindow;
%             currRegressionTable = cumStats{ind_subject}(ind_dof).regression_cumulative_singleWindow;
%             
%             [currFeatureTableSeg, currFeatureTableRest, segMask, restMask] = sepSegRest(currFeatureTable);
%             [currRegressionTableSeg, currRegressionTableRest, segMask, restMask] = sepSegRest(currIndivTable);
% %             [currRegressionTableSeg, currRegressionTableRest, segMask, restMask] = sepSegRest(currRegressionTable);
%             
%             % pull out the Rsq2
%             for ind_features = 1:length(allFeaturesSingle)
%                 currFeature = allFeaturesSingle{ind_features};
%                 
%                 regressionTableInd = find(strcmpi(currRegressionTableSeg.statType, currFeature));
%                 regressionTableSub = currRegressionTableSeg(regressionTableInd, :);
%                 featureTableSegTimeSingle{ind_dof, ind_features, ind_subject} = currFeatureTableSeg.time(2:end);
%                 featureTableSegDataSingle{ind_dof, ind_features, ind_subject} = regressionTableSub.Rsq2;
%                 featureTableSegMeanSingle{ind_dof, ind_features, ind_subject} = currFeatureTableSeg.mean;
%                 
% %                 if mean(currFeatureTableSeg.mean) < mean_mag_threshold 
% %                      featureTableSegDataSingle{ind_dof, ind_features, ind_subject} = 0*ones(size(featureTableSegDataSingle{ind_dof, ind_features, ind_subject}));
% %                 end
%                 
%                 regressionTableInd = find(strcmpi(currRegressionTableRest.statType, currFeature));
%                 regressionTableSub = currRegressionTableRest(regressionTableInd, :);
%                 featureTableRestTimeSingle{ind_dof, ind_features, ind_subject} = currFeatureTableRest.time;
%                 featureTableRestDataSingle{ind_dof, ind_features, ind_subject} = regressionTableSub.Rsq2;
%                 featureTableRestMeanSingle{ind_dof, ind_features, ind_subject} = currFeatureTableRest.mean;
%                 
% %                 if mean(currFeatureTableRest.mean) < mean_mag_threshold
% %                     featureTableRestDataSingle{ind_dof, ind_features, ind_subject} = 0*ones(size(featureTableSegDataSingle{ind_dof, ind_features, ind_subject}));
% %                 end
%             end
%         end
%     end
%     
%     plotStuff2('SingleSeg', allDofs, allFeaturesSingle, nSubject, featureTableSegTimeSingle, featureTableSegDataSingle, featureTableSegMeanSingle, outRegCsv, outputPathReg);
%     plotStuff2('SingleRest', allDofs, allFeaturesSingle, nSubject, featureTableRestTimeSingle, featureTableRestDataSingle, featureTableRestMeanSingle, outRegCsv, outputPathReg);
    
    
    % now assemble everything into a 3D array
    for ind_subject = 1:length(cumStats)
        currSubj = cumStats{ind_subject};
        
        if isempty(currSubj)
            continue
        end
        
        for ind_dof = 1:length(allDofs)
            if ind_dof > length(allDofs) / 2
                ind_dof_normUse = ind_dof - length(allDofs) / 2;
            else
                ind_dof_normUse = ind_dof;
            end
            
%             currStats = cumStats{ind_subject}(ind_dof);
            currFeatureTable = cumStats{ind_subject}(ind_dof).segStats_SingleWindow;
            currRegIndTable = cumStats{ind_subject}(ind_dof).regression_individual_singleWindow;
            currRegCumTable = cumStats{ind_subject}(ind_dof).regression_cumulative_singleWindow;
            
            currFeatureTable2 = cumStats{ind_subject}(ind_dof_normUse).segStats_SingleWindow;
            
            [currFeatureTableSeg, currFeatureTableRest, segMask, restMask] = sepSegRest(currFeatureTable);
            [currRegIndTableSeg, currRegIndTableRest, segMask, restMask] = sepSegRest(currRegIndTable);
            [currRegCumTableSeg, currRegCumTableRest, segMask, restMask] = sepSegRest(currRegCumTable);
            [currFeatureTableSeg2, currFeatureTableRest2, segMask, restMask] = sepSegRest(currFeatureTable2);
            
            for ind_features = 1:length(allFeaturesSingle)
                currFeature = allFeaturesSingle{ind_features};
                
                regressionTableSub = currRegIndTableSeg(strcmpi(currRegIndTableSeg.statType, currFeature), :);
                featureTableSegTimeSingle{ind_dof, ind_features, ind_subject} = currFeatureTableSeg.time;
                featureTableSegDataSingle{ind_dof, ind_features, ind_subject} = currFeatureTableSeg.(currFeature);
                featureTableSegRsq2Single{ind_dof, ind_features, ind_subject} = [0; regressionTableSub.Rsq2];
                featureTableSegMeanSingle{ind_dof, ind_features, ind_subject} = currFeatureTableSeg2.mean;

                regressionTableSub = currRegIndTableRest(strcmpi(currRegIndTableRest.statType, currFeature), :);
                featureTableRestTimeSingle{ind_dof, ind_features, ind_subject} = currFeatureTableRest.time;
                featureTableRestDataSingle{ind_dof, ind_features, ind_subject} = currFeatureTableRest.(currFeature);
                featureTableRestRsq2Single{ind_dof, ind_features, ind_subject} = [0; regressionTableSub.Rsq2];
                featureTableRestMeanSingle{ind_dof, ind_features, ind_subject} = currFeatureTableRest2.mean;
            end
        end
    end
    
%     allFeaturesMultiple = cumStats{1}(1).segStats_MultipleWindow.Properties.VariableNames(8:end);
%     % now assemble everything into a 3D array
%     for ind_subject = 1:length(cumStats)
%         currSubj = cumStats{ind_subject};
%         
%         if isempty(currSubj)
%             continue
%         end
%         
%         for ind_dof = 1:length(allDofs)
%             %             currStats = cumStats{ind_subject}(ind_dof);
%             currFeatureTable = cumStats{ind_subject}(ind_dof).segStats_MultipleWindow;
%             currRegressionTable = cumStats{ind_subject}(ind_dof).segStats_MultipleWindow;
%             [currFeatureTableSeg, currFeatureTableRest, segMask, restMask] = sepSegRest(currFeatureTable);
%             [currRegressionTableSeg, currRegressionTableRest, segMask, restMask] = sepSegRest(currRegressionTable);
%             
%             for ind_features = 1:length(allFeaturesMultiple)
%                 featureTableSegTimeMultiple{ind_dof, ind_features, ind_subject} = currFeatureTableSeg.time2;
%                 featureTableSegDataMultiple{ind_dof, ind_features, ind_subject} = currFeatureTableSeg.(allFeaturesMultiple{ind_features});
%                 
%                 featureTableRestTimeMultiple{ind_dof, ind_features, ind_subject} = currFeatureTableRest.time2;
%                 featureTableRestDataMultiple{ind_dof, ind_features, ind_subject} = currFeatureTableRest.(allFeaturesMultiple{ind_features});
%             end
%         end
%     end
    
    plotStuff('SingleSeg', allDofs, allFeaturesSingle, nSubject, featureTableSegTimeSingle, featureTableSegDataSingle, featureTableSegRsq2Single, featureTableSegMeanSingle, outStdCsv, outputPathStd);
    plotStuff('SingleRest', allDofs, allFeaturesSingle, nSubject, featureTableRestTimeSingle, featureTableRestDataSingle, featureTableRestRsq2Single, featureTableRestMeanSingle, outStdCsv, outputPathStd);
%     plotStuff('MultipleSeg', allDofs, allFeaturesMultiple, nSubject, featureTableSegTimeMultiple, featureTableSegDataMultiple, outCsv, outputPath);
%     plotStuff('MultipleRest', allDofs, allFeaturesMultiple, nSubject, featureTableRestTimeMultiple, featureTableRestDataMultiple, outCsv, outputPath);
end

function plotStuff(typeLabel, allDofs, allFeaturesSingle, nSubject, featureTableSegTime, featureTableSegData, featureTableRsq2, featureTableMean, outCsv, outputPath) 
    checkMkdir(outputPath);
    
    subjectColours = distinguishable_colors(nSubject);
    currInd = 0;
    
    % then plot everything
    for ind_dof = 1:length(allDofs)
        for ind_features = 1:length(allFeaturesSingle)
            currInd = currInd + 1;
            fprintf('Currently on %s - %u of %u\n', typeLabel, currInd, length(allDofs)*length(allFeaturesSingle));
            
            currDof = allDofs{ind_dof};
            currFeature = allFeaturesSingle{ind_features};
            figName = ['combined_' typeLabel '_' currDof '_' currFeature];
            figSavePath = fullfile(outputPath, figName);
            
            figName2 = ['indiv_' typeLabel '_' currDof '_' currFeature];
            figSavePath2 = fullfile(outputPath, figName2);
            
            figName3 = ['r2_' typeLabel '_' currDof '_' currFeature];
            figSavePath3 = fullfile(outputPath, figName3);
            
            h1 = figure('Position', [-1919 69 1920 964.8000]);
            h2 = figure('Position', [-1919 69 1920 964.8000]);
            h3 = figure('Position', [-1919 69 1920 964.8000]);
            hold on; 
            
            bSign = zeros(3, 1);
            
            for ind_subjects = 1:nSubject
                currColour = subjectColours(ind_subjects, :);
                currTime = featureTableSegTime{ind_dof, ind_features, ind_subjects};
                
                if isempty(currTime)
                    figure(h2);
                    ax1(ind_subjects) = subplot(3, 5, ind_subjects);
                    
                    figure(h3);
                    ax2(ind_subjects) = subplot(3, 5, ind_subjects);
                    continue
                end
                
                currData = featureTableSegData{ind_dof, ind_features, ind_subjects};
                currRsq2 = featureTableRsq2{ind_dof, ind_features, ind_subjects};
                currMean = featureTableMean{ind_dof, ind_features, ind_subjects};
                
                meanBlock(ind_subjects) = mean(currMean);
                
                if meanBlock(ind_subjects) < 1e-3
%                     lastR(ind_subjects) = 0;
%                     lastErr(ind_subjects) = 0;
                    
                    figure(h2);
                    ax1(ind_subjects) = subplot(3, 5, ind_subjects);
                    
                    figure(h3);
                    ax2(ind_subjects) = subplot(3, 5, ind_subjects);
                    continue
                end
                
                lastR(ind_subjects) = currRsq2(end);
                lastErr(ind_subjects) = currRsq2(end) - currRsq2(end-1);
                
                
                
                [b(:, ind_subjects), Rsq2(ind_subjects), X, yCalc2] = linearFit(currTime, currData);
                dataLen(ind_subjects) = length(currData);
                
                signageB = sign(b(2, ind_subjects));
                if signageB > 0
                    bSign(1) = bSign(1) + 1;
                elseif signageB < 0
                    bSign(2) = bSign(2) + 1;
                elseif isnan(b(2, ind_subjects))
                    bSign(3) = bSign(3) + 1;
                end
                
                currLabel = ['S' num2str(ind_subjects) ', b=' num2str(b(2, ind_subjects), '%0.4f'), ', R2=', num2str(Rsq2(ind_subjects), '%0.4f')];
                currLabel2 = ['S' num2str(ind_subjects), ...
                    ', R2=', num2str(lastR(ind_subjects), '%0.4f') ', R2E=', num2str(lastErr(ind_subjects), '%0.4f')];
                
                figure(h1);
                hold on;
                plot(currTime, currData, 'DisplayName', currLabel, 'Color', currColour, 'MarkerSize', 14, 'Marker', 'o', 'LineStyle', 'none');
                ph = plot(currTime, yCalc2, 'Color', currColour);
                ph.Annotation.LegendInformation.IconDisplayStyle = 'off';
                
                figure(h2);
                ax1(ind_subjects) = subplot(3, 5, ind_subjects);
                hold on;
                plot(currTime, currData, 'DisplayName', currLabel, 'Color', currColour, 'MarkerSize', 14, 'Marker', 'o', 'LineStyle', 'none');
                ph = plot(currTime, yCalc2, 'Color', [1 0 0]);
                ph.Annotation.LegendInformation.IconDisplayStyle = 'off';
                ylabel('Weight/dweight');
                title(currLabel);
                
                figure(h3);
                ax2(ind_subjects) = subplot(3, 5, ind_subjects);
                hold on;
                plot(currTime, currRsq2, 'DisplayName', currLabel, 'Color', currColour, 'MarkerSize', 14, 'Marker', 'o', 'LineStyle', '-');
%                 ph = plot(currTime, yCalc2, 'Color', [1 0 0]);
%                 ph.Annotation.LegendInformation.IconDisplayStyle = 'off';
                title(currLabel2);
                ylabel('Rsq');
%                 xlim('');
                ylim([0 1]);
            end
            
            linkaxes(ax1, 'y');
            linkaxes(ax2, 'y');
            
            % what is the mean/std slope and Rsq value? 
            meanB = mean(b(2, :));
            stdB = std(b(2, :));
            meanRsq = mean(Rsq2);
            stdRsq = std(Rsq2);
            
            figure(h1);
            titleStr = [typeLabel '_' currDof '_' currFeature, ...
                ', b=' num2str(meanB, '%0.4f') '\pm' num2str(stdB, '%0.4f'), ...
                ', bsign=(+)' num2str(bSign(1), '%0.4f') ', (-)' num2str(bSign(2), '%0.4f'), ...
                ', R2=', num2str(meanRsq, '%0.4f') '\pm' num2str(stdRsq, '%0.2f')];
            title(titleStr);
            
            legend('show');
            
            saveas(h1, figSavePath, 'png');
            saveas(h1, figSavePath, 'fig');
            close(h1);
            
            saveas(h2, figSavePath2, 'png');
            saveas(h2, figSavePath2, 'fig');
            close(h2);
            
            saveas(h3, figSavePath3, 'png');
            saveas(h3, figSavePath3, 'fig');
            close(h3);
            
            if ~exist(outCsv, 'file')
                header = 'typeLabel,dof,feature,b_mean,b_std,R2_mean,R2_std,bsign_plus,bsign_minus,bsign_nan';
                for ind_subjects = 1:nSubject
                    header = [header ',bmag_s' num2str(ind_subjects)];
                end
                for ind_subjects = 1:nSubject
                    header = [header ',r2_s' num2str(ind_subjects)];
                end
                for ind_subjects = 1:nSubject
                    header = [header ',r2e_s' num2str(ind_subjects)];
                end
                for ind_subjects = 1:nSubject
                    header = [header ',len_s' num2str(ind_subjects)];
                end
                for ind_subjects = 1:nSubject
                    header = [header ',mean_s' num2str(ind_subjects)];
                end
                header = [header '\n'];
            else
                header = '';
            end
            
            fid = fopen(outCsv, 'a');
            fprintf(fid, [header]);
            
            fprintf(fid, '%s,%s,%s', typeLabel, currDof, currFeature);
            fprintf(fid, ',%f,%f', meanB, stdB);
            fprintf(fid, ',%f,%f', meanRsq, stdRsq);
            fprintf(fid, ',%f,%f,%f', bSign(1), bSign(2), bSign(3));
            
            for ind_subjects = 1:nSubject
                fprintf(fid, ',%f', b(2, ind_subjects));
            end
            for ind_subjects = 1:nSubject
                fprintf(fid, ',%f', lastR(ind_subjects));
            end
            for ind_subjects = 1:nSubject
                fprintf(fid, ',%f', lastErr(ind_subjects));
            end
            for ind_subjects = 1:nSubject
                fprintf(fid, ',%f', dataLen(ind_subjects));
            end
            for ind_subjects = 1:nSubject
                fprintf(fid, ',%f', meanBlock(ind_subjects));
            end
            fprintf(fid, '\n');
            fclose(fid);
        end
    end    
end