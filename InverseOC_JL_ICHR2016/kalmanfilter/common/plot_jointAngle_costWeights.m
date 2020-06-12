function h = plot_jointAngle_costWeights(outputData, optimalControlInstance, featureInstance)
    h = figure;
    
    jointAngleNames = {featureInstance.model.joints.name};
    costFunctionNames = {optimalControlInstance.costFunctionStruct.name};
    
    [time, q, weight_averaged] = meanJointAngles(outputData);
    
    ax(1) = subplotJointAngles(subplot(311), time, q, 'Input Trajectory', jointAngleNames);
    plotBoxes(featureInstance.segments);
    
    time = zeros(size(outputData));
    weight_nonaveraged = zeros(numel(time), size(outputData(1).c, 2));
    resnormArray = zeros(size(outputData));
    timePts(1) = floor(length(outputData(1).time)/2) - floor(optimalControlInstance.windowAdvance/2);
    timePts(2) = floor(length(outputData(1).time)/2) + floor(optimalControlInstance.windowAdvance/2);
    
    for i = 1:length(outputData)
        ind = ((i-1)*2+1):((i)*2);
        
        time(ind) = outputData(i).time(timePts);
        weight_nonaveraged(ind, :) = [outputData(i).c; outputData(i).c];
        resnormArray(ind) = [outputData(i).resnorm outputData(i).resnorm];
    end
    
    ax(2) = subplotCostWeights(subplot(312), time, weight_nonaveraged, 'IOC', costFunctionNames);
%     ax(2) = subplotCostWeights(subplot(312), time, weight_averaged, 'IOC', costFunctionNames);
    plotBoxes(featureInstance.segments, [1 0 0], 0, -0.1, 1.1);
    ylim([-0.2 1.2]);
    
    ax(3) = subplotResnorm(subplot(313), time, resnormArray, 'Residual Norm');
%     plotBoxes(featureInstance.segments);
     
    linkaxes(ax, 'x');
end

function ax = subplotJointAngles(subplotHandle, time, data, plotTitle, legendNames)
    ax = subplot(subplotHandle);
    plot(time, data);
    title(plotTitle);
    ylabel('Joint angle [rad]');
%     xlabel('Time [s]');
    legend(legendNames);
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

function ax = subplotResnorm(subplotHandle, time, data, plotTitle)
    ax = subplot(subplotHandle);
    plot(time, data);
    title(plotTitle);
    ylabel('Resnorm');
    xlabel('Time [s]');
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