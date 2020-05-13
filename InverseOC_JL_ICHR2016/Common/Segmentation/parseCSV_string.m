function dataOut = parseCSV_segGuidance(path, settings)

if ~exist('settings', 'var')
    settings.skipLastEntry = 1;
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

numerical_index = [];
string_index = [];
for i=1:numel(headers)
    if ~isempty(strfind(headers{i},'Name'))
        string_index(end+1) = i;
    else
        numerical_index(end+1) = i;
    end
end

num_els = numel(headers);
data_num = zeros(numel(flines)-1,length(numerical_index));
data_str = cell(numel(flines)-1,length(string_index));

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
    
    strSplitArray = strsplit(flines{i},',');
    
    for j = 1:length(numerical_index)
        data_num(i-1,j) = str2double(strSplitArray{numerical_index(j)});
    end
    for j = 1:length(string_index)
        data_str{i-1,j} = strSplitArray{string_index(j)};
    end
    
    row = row+1;
end
data_num = data_num(1:row,:);
cells = mat2cell(data_num,row,ones(1,numel(headers(numerical_index))));
dataOut = cell2struct(cells,headers(numerical_index),2);

for i = 1:size(string_index,1)
    dataOut.(headers{string_index}) = data_str(1:row, i);
end

fclose(fid);
end

