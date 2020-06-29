function returnTime = parseTime(timeStamp)
    % given a flat string format (ie Cortex SDK), returns the time in seconds
    hrIndStart = 1;
    minIndStart = hrIndStart + 3;
    secIndStart = minIndStart + 3;
    
    hour = str2num(timeStamp(hrIndStart:hrIndStart+1));
    minutes = str2num(timeStamp(minIndStart:minIndStart+1));
    seconds = str2num(timeStamp(secIndStart:end));
    
    returnTime = hour*60*60 + minutes*60 + seconds; % convert to seconds
end