function writeCSV(file,data)
    fid = fopen(file,'w');

    names = fieldnames(data);
    r = numel(data.(names{1}));
    
    uncalibrated_index = zeros(1,numel(names));
    for j=1:numel(names)
        if(~isempty(strfind(names{j},'Uncalibrated')))
            data.(names{j}) = typecast(int32(data.(names{j})),'uint32');
            uncalibrated_index(j) = 1;
        end
    end
    
    c = struct2cell(data);
    mu = cell2mat(c(uncalibrated_index==1)');
    m = cell2mat(c(uncalibrated_index==0)');
    
    format = cell(1,numel(names));
    if(isempty(mu))
        format(:) = {'%f,'};
    elseif(isempty(m))
        format(:) = {'0x%08x,'};
    else
        format(1:size(m,2)) = {'%f,'};
        format(size(m,2)+1:size(m,2)+size(mu,2)) = {'0x%08x,'};
    end
        format = [format{:} '\r\n'];
    
    fprintf(fid,'%s%s\r\n',sprintf('%s,',names{uncalibrated_index==0}),sprintf('%s,',names{uncalibrated_index==1}));
    if(isempty(mu))
        for i=1:r
            fprintf(fid,format,m(i,:));
        end
    elseif(isempty(m))
        for i=1:r
            fprintf(fid,format,mu(i,:));
        end
    else
        for i=1:r
            fprintf(fid,format,m(i,:),mu(i,:));
        end
    end
    
    fclose(fid);
end
