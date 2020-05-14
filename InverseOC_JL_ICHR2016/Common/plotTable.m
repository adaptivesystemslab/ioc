function plotTable( arrayToPlot )
    % PURPOSE - takes in arrayToPlot and plots it
    % assumes first column is the x-axis (so time or frame rate)
    
    arrayLen = size(arrayToPlot, 2) - 1;
    
    figHan = figure;
    hold on;
    
    for x = 1:arrayLen
        plot(arrayToPlot(:,1), arrayToPlot(:,x));
        hold all;
    end
    
end