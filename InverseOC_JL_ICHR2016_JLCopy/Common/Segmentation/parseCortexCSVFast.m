function data = parseCortexCSVFast(path,start_line)

    if nargin < 2
       start_line = 2; 
    end

    fid = fopen(path);
    fdata = fread(fid,inf,'uint8=>char');
    flines = regexp(fdata','\r\n','split');

    headers = regexp(lower(flines{1}),',','split');

    end_comma = false;
    if(flines{1}(end) == ',')
        end_comma = true;
        headers = headers(1:end-1);    
    end
    
    %Trim
    headers = cellfun(@strtrim,headers,'UniformOutput', false);
    for i=1:numel(headers)
       headers{i} = strrep(headers{i},'.','_'); 
    end
    
    format = [];
    timestamp_index = [];
    for i=1:numel(headers)
        if ~isempty(strfind(headers{i},'timestamp'))
            format = strcat(format,'%12c,');
            timestamp_index = i;
        else
            format = strcat(format,'%f,');
        end
    end
    
    if ~end_comma
        format = format(1:end-1);
    end

    flines = flines(start_line:end);
    
    num_els = numel(headers);
    data_mat = zeros(numel(flines),num_els);

    for i=1:numel(flines)
        row = sscanf(flines{i},format);
        %Convert Timestamp to Unix
        time = datenum(char(row(timestamp_index:timestamp_index+11))','HH:MM:SS.FFF');
        data_mat(i,[1:timestamp_index-1 timestamp_index+1:end]) = row([1:timestamp_index-1 timestamp_index+12:end]);
        data_mat(i,timestamp_index) = time;
    end
    cells = mat2cell(data_mat,numel(flines),ones(1,numel(headers)));
    data = cell2struct(cells,headers,2);
    fclose(fid);

end