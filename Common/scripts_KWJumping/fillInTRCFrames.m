function trcOut = fillInTRCFrames(trc)

markerNames = fieldnames(trc.data);
markerNames = markerNames(3:end); % first 2 names are Frame # and Time

% If markers missing at beginning/end, copy first/last positions into missing frames
for m = 1:length(markerNames)
    % fill in initial frames
    if(trc.data.(markerNames{m})(1,:) == [0,0,0]) 
        idx_start = find(trc.data.(markerNames{m})(:,1) ~= 0,1);
        trc.data.(markerNames{m})(1:(idx_start-1),:) = repmat(trc.data.(markerNames{m})(idx_start,:),(idx_start-1),1);
    end
    % fill in final frames
    if(trc.data.(markerNames{m})(end,:) == [0,0,0])
        idx_end = find(trc.data.(markerNames{m})(:,1) ~= 0,1,'last');
        trc.data.(markerNames{m})((idx_end+1):end,:) = repmat(trc.data.(markerNames{m})(idx_end,:),(trc.OrigNumFrames - idx_end),1);
    end  
end

trcOut = trc;