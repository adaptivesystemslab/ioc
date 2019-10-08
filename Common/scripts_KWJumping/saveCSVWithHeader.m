function saveCSVWithHeader(filePath,header,data)

if(~contains(filePath,'.csv'))
    filePath = strcat(filePath,'.csv');
end

% write header
if(~isempty(header))
    if(isequal(class(header),'cell'))
        headerLast = header(end);
        header = [header(1:end-1);repmat({','},1,numel(header)-1)]; %insert commas
        header = header(:)';
        header = cell2mat(header);
        header = strcat(header,headerLast); % so there is no delimeter at the end of the header line
    end
    
    fid = fopen(filePath,'w');
    fprintf(fid,'%s\n',header);
    fclose(fid);
end

% add data
dlmwrite(filePath,data,'-append');