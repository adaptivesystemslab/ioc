% dataTRO plotting script
motionsofinterest = 1:2;

currTimeWindow = localSegTimeStart(motionsofinterest(1)):localSegTimeEnd(motionsofinterest(end));
sigDof = sigDofStdDev(subjQ(currTimeWindow, :), 3);
xlimWin = [currTimeWindow(1) - 1 currTimeWindow(end) + 1];

figure(h1);
clf
subplot(2, 1, 1)
plot(rawTime, subjQ(:, :));
plotBoxes(h1, subjSegTime.startTimeVal, subjSegTime.endTimeVal);

subplot(2, 1, 2)
plot(rawTime, subjdQ(:, :));
plotBoxes(h1, subjSegTime.startTimeVal, subjSegTime.endTimeVal);

figure(h2);
clf
for inda = 1:length(sigDof)
    switch length(sigDof)
        case 1
            % nothing
        case 2
            subplot(2, 1, inda);
        case 3
            subplot(3, 1, inda);
        case 4
            subplot(2, 2, inda);
        case 5
            subplot(2, 3, inda);
        case 6
            subplot(2, 3, inda);
    end
    
    activeSigDof = sigDof(inda);
    
    hold on
    
    skimmedDownInd = 1:10:length(subjQ(:, 1));
    scatter(subjQ(skimmedDownInd, activeSigDof), subjdQ(skimmedDownInd, activeSigDof), 'b.'); hold on
    
    title([num2str(sigDof(inda)) ': ' TRO2009_DOFName(sigDof(inda))]);
    
    
end

for timeWindowInd = motionsofinterest
    currTimeWindow = localSegTimeStart(timeWindowInd):localSegTimeEnd(timeWindowInd);
    figure(h1);
    
    switch timeWindowInd
        case 1
            colour = 'ro';
        case 2
            colour = 'go';
        case 18
            colour = 'ko';
        case 19
            colour = 'mo';
    end
    
    subplot(2, 1, 1)
    hold on
    plot(rawTime(currTimeWindow), subjQ(currTimeWindow, sigDof), colour)
    xlim(xlimWin)
    
    subplot(2, 1, 2)
    hold on
    plot(rawTime(currTimeWindow), subjdQ(currTimeWindow, sigDof), colour)
    xlim([38 41])
    
    ylim([-5 5]);
    for indblah = 1:length(subjSegTime.segmentName)
        if subjSegTime.startTimeVal(indblah) > xlimWin(1) && ...
                subjSegTime.startTimeVal(indblah) < xlimWin(2)
            if mod(indblah,2)
                ypos = -4;
            else
                ypos = -3;
            end
            text(subjSegTime.startTimeVal(indblah), ypos, subjSegTime.segmentName(indblah));
        end
    end
    hold off
    
    figure(h2);
    for inda = 1:length(sigDof)
        switch length(sigDof)
            case 1
                % nothing
            case 2
                subplot(2, 1, inda);
            case 3
                subplot(3, 1, inda);
            case 4
                subplot(2, 2, inda);
            case 5
                subplot(2, 3, inda);
            case 6
                subplot(2, 3, inda);
        end
        
        activeSigDof = sigDof(inda);
        
        firstHalfWindow = currTimeWindow(1:ceil(length(currTimeWindow)/2));
        secondHalfWindow = currTimeWindow(floor(length(currTimeWindow)/2):length(currTimeWindow));
        
        hold on
        
        scatter(subjQ(currTimeWindow(end-3:end), activeSigDof), subjdQ(currTimeWindow(end-3:end), activeSigDof), [colour(1) 'x']);
        scatter(subjQ(secondHalfWindow, activeSigDof), subjdQ(secondHalfWindow, activeSigDof), colour);
        scatter(subjQ(currTimeWindow(1:3), activeSigDof), subjdQ(currTimeWindow(1:3), activeSigDof), [colour(1) 'x']);
        scatter(subjQ(firstHalfWindow, activeSigDof), subjdQ(firstHalfWindow, activeSigDof), colour);
        
    end
end