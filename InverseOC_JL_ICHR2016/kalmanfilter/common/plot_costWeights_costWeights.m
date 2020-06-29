function h = plot_costWeights_costWeights(outputDataDoc, outputDataIoc, costWeightNames)
    h = figure;
    
    [timeDoc, ~, weightDoc] = meanJointAngles(outputDataDoc);
    [timeIoc, ~, weightIoc] = meanJointAngles(outputDataIoc);
    
    ax(1) = subplotCostWeights(subplot(211), timeDoc, weightDoc, 'DOC', costWeightNames);
    ax(2) = subplotCostWeights(subplot(212), timeIoc, weightIoc, 'IOC', costWeightNames);
    
    linkaxes(ax);
end

function ax = subplotCostWeights(subplotHandle, time, data, plotTitle, legendNames)
    ax = subplot(subplotHandle);
    plotVal = [];
    for i = 1:size(data, 1)
        plotVal(i, :) = data(i, :)/sum(data(i, :));
    end
    area(time, plotVal);
    title(plotTitle);
    ylabel('Weights normalized');
%     xlabel('Time [s]');
%     ylim([-0.2 1.2]);
    legend(legendNames);
end

function [time, data, weight] = meanJointAngles(outputData)
    time = unique([outputData.time]);
    data = zeros(numel(time), size(outputData(1).q, 2));
    weight = zeros(numel(time), size(outputData(1).c, 2));
    n = zeros(size(time));
    
    for ind_outputData = 1:length(outputData)
        for ind_time = 1:length(time)
              findInd = find(time(ind_time) == outputData(ind_outputData).time);
              
              for ind_findInd = 1:length(findInd)
                  currInd = findInd(ind_findInd);
                  data(ind_time, :) = data(ind_time, :) + outputData(ind_outputData).q(currInd, :);
                  weight(ind_time, :) = weight(ind_time, :) + outputData(ind_outputData).c;
                  n(ind_time) = n(ind_time) + 1;
              end
        end
    end
    
    for ind_n = 1:length(n)
        data(ind_n, :) = data(ind_n, :) / n(ind_n);
        weight(ind_n, :) = weight(ind_n, :) / n(ind_n);
    end
end