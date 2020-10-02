function data = parseCortexCSV(file, shimmerParam, startCrop, endCrop)

    if ~exist('startCrop', 'var')
        startCrop = 1;
    end
    
    if ~exist('endCrop', 'var')
        endCrop = 1;
    end

    fid = fopen(file);
    %Read First Line
    tline = fgets(fid);

    %Save Headers
    % remove stray commas and other odd spacing at the end of the headers
    commaCheckLoop = 1;
    while commaCheckLoop
        if isempty(regexp(tline(end),'[A-Za-z0-9]','end'))
            tline = tline(1:end-1);
        else
            commaCheckLoop = 0;
        end
    end
    
    headers = regexp(tline,',','split');
    %Trim
    headers = cellfun(@strtrim,headers,'UniformOutput', false);
    for i=1:numel(headers)
        if strcmpi(headers{i}, 'Timestamp')
            headers{i} = 'SystemMSTimeStamp';
        end
        
       headers{i} = strrep(headers{i},'.','_'); % parseCSV doesn't like dots
       headers{i} = strrep(headers{i},'-','_'); % or dashes
       
       % check the contents of the headers, if empty
    end
    
    num_els = numel(headers);
    %count number of lines
    %Number of lines in the file excluding the header
    nLines = 0;
    while (fgets(fid) ~= -1),
      nLines = nLines+1;
    end
    
    cells = cell(1,num_els);
    for i=1:num_els
        cells{i} = zeros(nLines-startCrop-endCrop+1,1);
    end
    data = cell2struct(cells,headers,2); 
    
    %Move seek to second line
    fseek(fid,0,'bof');
    %Skip header string again 
    fgets(fid);
    for i = 1:nLines-endCrop % skip the last line because sometime the cortex recorder cuts out before the last line finishes
        tline = fgets(fid);
        
        if i < startCrop
            % skip the first 'startCrop' amount
            continue
        end
        
        str_data = regexp(tline,',','split');
        for j=1:num_els
            if(~isempty(strfind(str_data{j},':')))
               %Convert Timestamp to Unix
               m_date_collected = shimmerParam.date;
               
%                timeRawOrig = datenum(str_data{j},'HH:MM:SS.FFF')*86400;
               timeRaw = parseTime(str_data{j});
               
               timeNow = timeRaw + date2utc(m_date_collected,0); % TODO may need to shift by 4 hours for UTC
               cortexTime = timeNow * 1000;
               data.(headers{j})(i-startCrop+1) = floor(cortexTime);
%                data.(headers{j})(i) = floor(86400*1000 * (time - datenum('01-Jan-1970')));
            else
                %Data String
                data.(headers{j})(i-startCrop+1) = sscanf(str_data{j},'%f');
            end
        end
    end
    
    % close the handles
    fclose(fid);
end


%     m_date_collected = shimmerParam.date;
% %     ack = datenum(startDate,'yyyy_mm_dd_HH_MM_SS');
%     segTimeRaw = segmentInfo.t;
% 
%     % reparse it as 'time'
%     if ~iscell(segTimeRaw)
%         segTimeNow = segTimeRaw + date2utc(m_date_collected,0);
%         segTime = segTimeNow * 1000;
%     else
%         for i = 1:length(segTimeRaw)
% %                     time = datenum(str_data{j},'HH:MM:SS.FFF');
% %                old_data.(headers{j})(i) = floor(86400*1000 * (time - datenum('01-Jan-1970')));
%                
%             segTimeNow = segTimeRaw{i} + date2utc(m_date_collected,0); % time since 1970 in seconds
%             segTime{i} = segTimeNow * 1000; % time since 1970 in ms
%             
% %             blah = datestr(segTime{i}/(60*60*24*1000) + datenum('01-Jan-1970'))
%         end
