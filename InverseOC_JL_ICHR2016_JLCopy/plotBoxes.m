function [minY, maxY, h_line] = plotBoxes(h, currTimeStart, currTimeEnd, colorToUse, offset, maxY, minY)
    % plotBoxes(h, currTimeStart, currTimeEnd, colorToUse, offset, maxY, minY)
    % plots the segmentation windows on an existing figure handle
    
%     figure(h);
  
    ylimVal = ylim;
    
    if ~exist('minY', 'var')
%         minY = ylimVal(1) * 1.0005;
%         maxY = ylimVal(2)* 1.0005;
        maxYLim = max(abs(ylimVal))*1.0005;
        
        minY = -maxYLim;
        maxY =  maxYLim;
    end
    
    if ~exist('colorToUse', 'var')
        colorToUse{1} = [1 0 0];
    elseif ~iscell(colorToUse)
        colorToUseA{1} = colorToUse;
        colorToUse = colorToUseA;
    end
    
    if ~exist('offset', 'var')
        offset = 0;
    end
    
    hold on

    for i = 1:length(currTimeStart)
%         colorToUse = 'r'; % orig
        % flip around the offset so it's easier to see overlaping boxes
        offset = offset * -1;
        colorCounter = colorToUse{mod(i, length(colorToUse)) + 1};
        
        startMarker = currTimeStart(i);
        endMarker = currTimeEnd(i);
        
        if i == 1
            h_line = plot([startMarker startMarker], [minY maxY] + offset, 'Color', colorCounter);
        else
            plot([startMarker startMarker], [minY maxY] + offset, 'Color', colorCounter);
        end
        
        plot([endMarker endMarker], [minY maxY] + offset,  'Color', colorCounter);
        plot([startMarker endMarker], [minY minY] + offset,  'Color', colorCounter);
        plot([startMarker endMarker], [maxY maxY] + offset,  'Color', colorCounter);
    end
    
%     ylim([minY maxY] * 1.1); 

    hold off
end