function al
    clearvars
    clc

    addpath(genpath(fullfile('Common')));
    addpath(genpath('Libraries/rl/ik_framework/common'));
    addpath(genpath('Libraries/rl/ik_framework/instance_expressive'));
    addpath(genpath('Libraries/rl/General_FKEKF/DynamicsModelMatlab/MatlabWrapper'));

     baseSourceFolder = '../results/ioc';
    %baseSourceFolder = 'D:/results/expressive_ioc/PamelasIOCResults/';
    commentsArray = {''};
    SVDArray = {'10000'};
    thresholdMultiplierArray = [1];

    nowStrMaster = datestr(now, 'yyyymmddHHMMSS');

    for ind_multiplierarray = 1:length(thresholdMultiplierArray)    
        for ind_commentsarray = 1:length(commentsArray)
            for ind_svdarray = 1:length(SVDArray)

                nowStr = datestr(now, 'yyyymmddHHMMSS');
                comments = commentsArray{ind_commentsarray};
                svdEntry = SVDArray{ind_svdarray};
                thresholdMultiplier = thresholdMultiplierArray(ind_multiplierarray);

                basePath_array =           {};
                basePathPlot_array =       {};
                basePathMasterPlot_array = {};        

                basePath_array{1} =           [];
                basePathPlot_array{1} =       [];
                basePathMasterPlot_array{1} = [];

                basePath_array{end+1} =           [baseSourceFolder 'expressive_ioc/'];
                basePathPlot_array{end+1} =       [baseSourceFolder 'expressive_ioc/' 'p_' nowStr '_' comments '_' svdEntry '/'];
                basePathMasterPlot_array{end+1} = basePathPlot_array{end};

                for i = 1:length(basePath_array)
                    if isempty(basePath_array{i})
                        continue
                    end

                    basePath = basePath_array{i};
                    basePathPlot = basePathPlot_array{i};
                    basePathMasterPlot = basePathMasterPlot_array{i};

                    masterSource = [basePathPlot 'log_overall.csv'];
                    masterFigureSave = [basePathPlot '/fig/'];

                    checkMkdir(basePathMasterPlot);
                    checkMkdir(masterFigureSave);

                    dir_layer1 = dir(basePath);

                    counter = 0;
                    totalIndWeighting = [];
                    accLess = {};
                    accMore = {};
                    maxIndTrack = [];
                    maxAccIndiv = [];
                    instanceNameCounter = {};
                    % rmse_reportCounter = [];
                    windowedRmse_mean = [];
                    windowedRmse_std = [];
                    windowedResnorm_mean = [];
                    windowedResnorm_std = [];

                    currDataset = strsplit(basePath_array{i}, '/'); 
                    masterSaveLog = fullfile(basePathMasterPlot, ['outputFile_plot_' currDataset{end-2} '_' nowStr '_' comments '_' svdEntry '.csv']);
                    masterSaveLogSeg = fullfile(basePathMasterPlot, ['outputFile_segment_' currDataset{end-2} '_' nowStr '_' comments '_' svdEntry '.csv']);

                    masterSaveOverall = fullfile(baseSourceFolder, ['overallFile_' nowStrMaster '.csv']);
                    masterSaveOverallSeg = fullfile(baseSourceFolder, ['overallFile_segment_' nowStrMaster '.csv']);

                    jumpInd = 0;

                    for ind_layer1 = 1:length(dir_layer1)
                        currDir_layer1 = dir_layer1(ind_layer1);
                        currDir_layer1_fullName = [basePath currDir_layer1.name '/'];

                        splitTheStr = strsplit(currDir_layer1.name, '_');

                        if ~currDir_layer1.isdir
                            continue
                        elseif length(currDir_layer1.name) < 3
                            continue
                        end

                        if length(splitTheStr) > 8 && ~strcmpi(splitTheStr{9}, svdEntry)
                            continue
                        end
                         
                        dir_layer2 = dir(currDir_layer1_fullName);
                        for ind_layer2 = 1:length(dir_layer2)
                            currDir_layer2 = dir_layer2(ind_layer2);
                            currDir_layer2_fullName = [currDir_layer1_fullName currDir_layer2.name '/'];

                            if strcmpi(currDir_layer2.name, 'log_overall.csv')
                                % copy it to the master source
                                fId1 = fopen([currDir_layer2_fullName(1:end-1)]);
                                fId2 = fopen(masterSource, 'a+');

                                tline = fgetl(fId1);
                                while ischar(tline)
                                    fprintf(fId2, '%s\n', tline);
                                    tline = fgetl(fId1);
                                end

                                fclose(fId1);
                                fclose(fId2);

                                continue
                            elseif ~currDir_layer2.isdir
                                continue
                            elseif length(currDir_layer2.name) < 3
                                continue
                            end

%                             if ~strcmp(currDir_layer2.name, '20190725_074549')
%                                 continue
%                             end
                            
                            dir_layer3 = dir([currDir_layer2_fullName]);

                            for ind_layer3 = 1:length(dir_layer3)
                                currDir_layer3 = dir_layer3(ind_layer3);
                                currDir_layer3_fullName = [currDir_layer2_fullName currDir_layer3.name '/'];

                                if ~currDir_layer3.isdir
                                    continue
                                elseif length(currDir_layer3.name) < 3
                                    continue
                                end

                                fprintf('Loading %s\n', currDir_layer3_fullName);

                                dir_layer4 = dir([currDir_layer3_fullName '*end.mat']);
                                matFilesToLoad = {};
                                indStart = [];

                                if isempty(dir_layer4)
                                    continue
                                end

                    % % %             try
                                 for ind_layer4 = 1:length(dir_layer4)
                                     currDir_layer4 = dir_layer4(ind_layer4);
                                     matFilesToLoad{ind_layer4} = [currDir_layer3_fullName currDir_layer4.name];

                                     subtextString = strsplit(currDir_layer4.name, '_');
                                     indStart(ind_layer4) = str2num(subtextString{3});
                                     indEnd(ind_layer4) = str2num(subtextString{4});
                                 end

                                 [Y, I] = sort(indStart);
                                 indStart = indStart(I);
                                 indEnd = indEnd(I);
                                 indTotal = [indStart' indEnd'];
                                 matFilesToLoadUse = matFilesToLoad(I);

                                 instanceName = currDir_layer3.name;
                                 currDir_layer4_plotFullName = [basePathPlot currDir_layer1.name '/' currDir_layer3.name];

                                 if ~isempty(matFilesToLoadUse)    
                                     jumpInd = jumpInd + 1;
                                     masterSave1 = fullfile(basePathMasterPlot, [currDir_layer1.name '_fig_' nowStr]);
                                     [t_obs{jumpInd}, avgWeightArray_ioc{jumpInd}, avgRatioArray_ioc{jumpInd}, ...
                                         cost_function_names_sorted, outputStruct(jumpInd)] = loadFileAndPlot_overlapMode(matFilesToLoadUse, instanceName, indTotal, currDir_layer4_plotFullName, masterSaveOverall, masterSave1);
                                     writeTrajectoryToFile(currDir_layer4_plotFullName,instanceName, cost_function_names_sorted, t_obs{jumpInd}, avgWeightArray_ioc{jumpInd}, avgRatioArray_ioc{jumpInd});
                                 end

                                 fclose all;
                            end
                        end
                    end
                end
            end
        end
    end
    
    outputPath = '../../expressiveiocData/results/expressive_ioc/';
    
    for i = 1:length(outputStruct)
        time = outputStruct(i).indices_forward;
        avgRatio = outputStruct(i).ratio_forward';
        avgWeight = outputStruct(i).weights_forward';
        writeTrajectoryToFile(outputPath, ['forward_' nowStr '_' outputStruct(i).currInstName], cost_function_names_sorted, time, avgWeight, avgRatio);
%         
%         time = outputStruct(i).indices_back;
%         avgRatio = outputStruct(i).weights_back';
%         avgWeight = zeros(size(avgRatio));
%         writeTrajectoryToFile(outputPath, ['back_' nowStr '_' outputStruct(i).currInstName], cost_function_names_sorted, time, avgWeight, avgRatio);
    end
%     
%     
%     arousalScore = [0.522628122
% 0.110159458
% -0.388556087
% 0.335734153
% -0.347126836
% -0.776532296
% -1.053881223
% -0.311248702
% -0.316066468
% -0.998674452
% -0.693621625
% -0.312834034
% -0.303949431
% 0.032663715
% 0.10988623
% -0.477045052
% -0.939540501
% -0.972105287
% -1.287652728
% -1.033672727];
% 
% pleasureScore = [0.672994937
% 0.777853215
% -0.291668873
% 0.177429849
% 0.037442583
% -0.719379497
% -1.009978889
% -0.228532727
% -0.31612139
% -0.867479969
% -0.658264581
% -0.18745559
% -0.183213222
% 0.221503082
% -0.130131269
% -0.463576735
% -0.850457804
% -0.91907227
% -1.047399302
% -0.620772248]'
% 
% group = [3
% 3
% 2
% 3
% 2
% 1
% 1
% 2
% 2
% 1
% 1
% 2
% 2
% 2
% 3
% 1
% 1
% 1
% 1
% 1
% ];
% 
%     
%     weightMean{1} = [];
%     weightMean{2} = [];
%     weightMean{3} = [];
%     stdMean{1} = [];
%     stdMean{2} = [];
%     stdMean{3} = [];
%     for j = 1:size(outputStruct(1).weights_forward, 2)
%         cumArray{1} = [];
%         cumArray{2} = [];
%         cumArray{3} = [];
%         for i = 1:jumpInd
%             if sum(outputStruct(i).weights_forward(:, j)) > 0
%                 cumArray{group(i)} = [cumArray{group(i)} outputStruct(i).weights_forward(:, j)];
%             else
%                 fprintf('Rejecting %u, %u, g%u\n', j, i, group(i));
%             end
%         end
%         
%         for i = 1:3
%             weightMean{i} = [weightMean{i}; mean(cumArray{i}, 2)'];
%             stdMean{i} = [stdMean{i}; std(cumArray{i}, 0, 2)'];
%         end
%     end
%     
%     figure;
%     t = 0:0.05:1;
%     cfToPlot = [2 3 7 8 9 15];
%     
%     groupToPlot = 1;
%     for i = 1:length(cfToPlot)
%         cfToPlotCurr = cfToPlot(i);
%         axToPlot = i;
%         ax(axToPlot) = subplot(3, 6, axToPlot); 
%         barwitherr(stdMean{groupToPlot}(:, cfToPlotCurr), t, weightMean{groupToPlot}(:, cfToPlotCurr));
%         title(cost_function_names_sorted{cfToPlotCurr});
%     end
%     
%     groupToPlot = 2;
%     for i = 1:length(cfToPlot)
%         cfToPlotCurr = cfToPlot(i);
%         axToPlot = i + length(cfToPlot);
%         ax(axToPlot) = subplot(3, 6, axToPlot);
%         barwitherr(stdMean{groupToPlot}(:, cfToPlotCurr), t, weightMean{groupToPlot}(:, cfToPlotCurr));
%         title(cost_function_names_sorted{cfToPlotCurr});
%     end
%     
%     groupToPlot = 3;
%     for i = 1:length(cfToPlot)
%         cfToPlotCurr = cfToPlot(i);
%         axToPlot = i + 2*length(cfToPlot);
%         ax(axToPlot) = subplot(3, 6, axToPlot);
%         barwitherr(stdMean{groupToPlot}(:, cfToPlotCurr), t, weightMean{groupToPlot}(:, cfToPlotCurr));
%         title(cost_function_names_sorted{cfToPlotCurr});
%     end
%     
%     linkaxes(ax);
% 
%     
%     figure;
%     t = 0:0.05:1;
%     for i = 1:15
%         ax(i) = subplot(3, 5, i);
% %         plot(weightMean(:, i));
%         barwitherr(stdMean{3}(:, i), t, weightMean{1}(:, i));
%         title(cost_function_names_sorted{i});
%     end
%     linkaxes(ax);
end