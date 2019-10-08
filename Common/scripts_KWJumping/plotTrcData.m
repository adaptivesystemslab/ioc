function plotTrcData(trc)

markerNames = fieldnames(trc.data);

% Below is now done in main script
% % Rotate markers so subject jumps in positive X direction
% for m = 3:length(markerNames)
%     trc.data.(markerNames{m}) = (rotz(-pi/2)*trc.data.(markerNames{m})')';
% end


for i = 1:length(trc.data.Time)
    figure(1); clf; hold on; grid on;
    for m = 3:10 % torso markers
        pos = trc.data.(markerNames{m})(i,:);
        plot3(pos(1),pos(2),pos(3),'ko');
    end
    for m = 11:14 % left arm markers
        pos = trc.data.(markerNames{m})(i,:);
        plot3(pos(1),pos(2),pos(3),'b*');
    end
    for m = 15:18 % right arm markers
        pos = trc.data.(markerNames{m})(i,:);
        plot3(pos(1),pos(2),pos(3),'r*');
    end
    for m = 19:25 % left leg markers
        pos = trc.data.(markerNames{m})(i,:);
        plot3(pos(1),pos(2),pos(3),'bo');
    end
    for m = 26:length(markerNames) % right leg markers
        pos = trc.data.(markerNames{m})(i,:);
        plot3(pos(1),pos(2),pos(3),'ro');
    end
    
    xlabel('X'); ylabel('Y');
    axis([-500 1500 -1000 1000 0 2000]);
    view(60,15);
    
    pause(0.001); % required or else figure won't update for viewer
end
    