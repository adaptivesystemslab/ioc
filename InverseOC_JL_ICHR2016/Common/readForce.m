function force_struct = readForce(filename)

%Read some Force file properties
fid = fopen(filename);
linenum = 2;

force_struct = struct();
%Read lines 2:4 to get some info about force plates
line = textscan(fid, '%s', 6, 'delimiter', '=', 'headerlines', 1);
line = line{1};

fields = line(1:2:end);
values = str2double(line(2:2:end));

for i=1:numel(fields)
    force_struct.(fields{i}) = values(i);
end

%Now read the header information
line = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', 1);
fields = textscan(line{1}{1},'%s','delimiter','\t');
fields = fields{1};
fields(ismember(fields,'')) = [];
%Remove non alphanumeric characters
fields = regexprep(fields,'[^a-zA-Z0-9 -]','');

fclose(fid);

%Read the actual data
M = dlmread(filename,'\t',5,0);

for col=1:numel(fields);
    force_struct.data.(fields{col}) = M(:,col);
end

end