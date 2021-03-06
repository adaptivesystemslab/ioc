function status=writeTrc(trc_struct,filename)
%Writes a trc strucutre into a trc motion capture file 

%Open file for writing
fid = fopen(filename,'w');
%Write PathFileType and pathfilename
fprintf(fid,'PathFileType\t4\t(X/Y/Z)\t%s\t\r\n',filename);

%Write header fields and their values
fprintf(fid,'%s\t',trc_struct.header_fields{:});
fprintf(fid,'\r\n');
for i=1:numel(trc_struct.header_fields)
    if(strcmp(class(trc_struct.(trc_struct.header_fields{i})),'uint32'))
       fprintf(fid,'%u\t',trc_struct.(trc_struct.header_fields{i}));         
    elseif(strcmp(class(trc_struct.(trc_struct.header_fields{i})),'double'))
       fprintf(fid,'%.2f\t',trc_struct.(trc_struct.header_fields{i})); 
    else
       fprintf(fid,'%s\t',trc_struct.(trc_struct.header_fields{i})); 
    end
end
fprintf(fid,'\r\n');

%Print Data Headers
fprintf(fid,'%s\t',trc_struct.data_fields{1:2});
fprintf(fid,'%s\t\t\t',trc_struct.data_fields{3:end});
fprintf(fid,'\r\n');
%Print X1 Y1 Z1....XN YN ZN
%First we skip the Frame# and Time
fprintf(fid,'\t\t');
for i=1:trc_struct.NumMarkers
   fprintf(fid,'X%d\tY%d\tZ%d\t',i,i,i); 
end
fprintf(fid,'\r\n\r\n');
fclose( fid );

%Write Tab Delimited Data
trc_mat = zeros(trc_struct.NumFrames,2+trc_struct.NumMarkers*3);
col = 1;
for i=1:numel(trc_struct.data_fields)
   if ~isempty(strfind(trc_struct.data_fields{i},'#'))
       field = trc_struct.data_fields{i}(1:end-1);
   else
       field = trc_struct.data_fields{i};
   end 
   data =  trc_struct.data.(field);
   trc_mat(:,col:col+size(data,2)-1)=data;
   col = col+size(data,2);
end

dlmwrite(filename,trc_mat,'delimiter','\t','precision','%g', 'newline', 'pc', '-append');
end