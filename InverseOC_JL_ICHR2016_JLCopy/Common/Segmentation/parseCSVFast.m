%better read
function data = parseCSVFast(path)

fid = fopen(path);
fdata = fread(fid,inf,'uint8=>char');
flines = regexp(fdata','\r\n','split');

headers = regexp(lower(flines{1}),',','split');

end_comma = false;
if(flines{1}(end) == ',')
    end_comma = true;
    headers = headers(1:end-1);    
end

format = [];
uncalibrated_index = [];
for i=1:numel(headers)
    if ~isempty(strfind(headers{i},'uncalibrated'))
        format = strcat(format,'%li,');
        uncalibrated_index(end+1) = i;
    elseif ~isempty(strfind(headers{i},'calibrated'))
        format = strcat(format,'%f,');
    end
end

if ~end_comma
    format = format(1:end-1);
end

num_els = numel(headers);
data_mat = zeros(numel(flines)-1,num_els);

for i=2:numel(flines)
    data_mat(i-1,:) = sscanf(flines{i},format);
    data_mat(i-1,uncalibrated_index) = typecast(uint32(data_mat(i-1,uncalibrated_index)),'int32');
end
cells = mat2cell(data_mat,numel(flines)-1,ones(1,numel(headers)));
data = cell2struct(cells,headers,2);
fclose(fid);
end

