%better read
function data = parseCSV(path, settingsInput)

if ~exist('settingsInput', 'var')
    settings.skipLastEntry = 1;
elseif isnumeric(settingsInput) % legacy support
    settings.skipLastEntry = settingsInput;
else
    settings = settingsInput;
end

fid = fopen(path);

fdata = fread(fid,inf,'uint8=>char');

% if this breaks, try using {\r\n|\n} for expression
flines = regexp(fdata','(\r)?\n','split');

headers = regexp(flines{1},',','split');

end_comma = false;
if(flines{1}(end) == ',')
    end_comma = true;
    headers = headers(1:end-1);    
end

format = [];
uncalibrated_index = [];
for i=1:numel(headers)
    if ~isempty(strfind(headers{i},'SystemTimeStamp'))
        headers{i} = 'SystemTimeStamp'; % sometimes the header would have random char at the start
        format = strcat(format,'%f,');
    
    elseif ~isempty(strfind(headers{i},'Uncalibrated'))
        format = strcat(format,'%li,');
        uncalibrated_index(end+1) = i;
    
    elseif ~isempty(strfind(headers{i},'Timestamp'))
        if i == 2
            % account for old cortex SDK format
            format = strcat(format,'%[1234567890.:],');
        else
            format = strcat(format,'%f,');
        end
        
    else
        format = strcat(format,'%f,');
    end
end

if ~end_comma
    format = format(1:end-1);
end

num_els = numel(headers);
data_mat = zeros(numel(flines)-1,num_els);

row = 0;
if settings.skipLastEntry
    rowsToUse = numel(flines)-1; % sometimes the last line would terminate prematurely, so we should skip it in the parsing
else
    rowsToUse = numel(flines);
end

for i=2:rowsToUse 
    if(isempty(flines{i}))
       break; 
    end
    data_mat(i-1,:) = sscanf(flines{i},format);
    data_mat(i-1,uncalibrated_index) = typecast(uint32(data_mat(i-1,uncalibrated_index)),'int32');
    row = row+1;
end
data_mat = data_mat(1:row,:);
cells = mat2cell(data_mat,row,ones(1,numel(headers)));
data = cell2struct(cells,headers,2);
fclose(fid);
end

