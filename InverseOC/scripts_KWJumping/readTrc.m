function trc_struct = readTrc(filename,interpolate,new_rate)
%Reads a TRC motion capture file into a structure

    %Read the data part of TRC
    %The first two columns are frame number and time stamp the rest is data
    M = dlmread(filename,'\t',6,0);
    
    %Read some TRC file properties
    fid = fopen(filename);
    linenum = 2;
    line = textscan(fid, '%s', 1, 'delimiter', '\n', 'headerlines', linenum-1);
    fields = textscan(line{1}{1},'%s','delimiter','\t');
    fields = fields{1};
    fields(ismember(fields,'')) = [];
    line = textscan(fid, '%s', 1, 'delimiter', '\n');
    values = textscan(line{1}{1},'%s','delimiter','\t');
    values = values{1};
    values(ismember(values,'')) = [];
    %Put this information into our structure
    trc_struct = struct();
    trc_struct.header_fields = fields;
    for i=1:numel(fields)
        [num, status] = str2num(values{i});
        if(status)
           trc_struct.(fields{i}) = num;
        else
           trc_struct.(fields{i}) = values{i};
        end
    end
    rate = trc_struct.DataRate;
    
    
    %Now read the data names in the 4rth line of TRC file
    line = textscan(fid, '%s', 1, 'delimiter', '\n');
    fclose(fid);
    fields = textscan(line{1}{1},'%s','delimiter','\t');
    fields = fields{1};
    fields(ismember(fields,'')) = [];
    trc_struct.h = fields;
    %The first field is always Frame# we have to get rid of # to make it a
    %valid struct property
    fields{1} = fields{1}(1:end-1);
    %Sometimes we read a column of zeros at the end, get rid of it 
    M = M(:,1:(numel(fields)-2)*3+2);
    
    %Linearly Interpolate the missing markers
    %Missing markers are 0 due to dlmread
    if nargin > 1 && interpolate == 1
        for col=3:(numel(fields)-2)*3+2;
           marker_data = M(:,col);
           times = M(marker_data ~= 0,2);
           if sum(marker_data ~= 0) > 10
              M(:,col) = interp1(times,marker_data(marker_data ~= 0),M(:,2)); 
           else
              disp('asdf'); 
           end
        end
    end
    
    %Resample if the sampling rate was incorrect
    if nargin > 2 && new_rate ~= trc_struct.DataRate
        new_times = linspace(0,M(end,2),(M(end,2)+2/new_rate)/(1/new_rate));
        for col=3:(numel(fields)-2)*3+2;
           marker_data = M(:,col);
           times = M(:,2);
           M_new(:,col-2) = interp1(times,marker_data(marker_data ~= 0),new_times); 
        end
        
        M_new = [(1:size(M_new,1))' new_times' M_new];
        M = M_new;
        rate=new_rate;
    end
    
    
    %Now we assign the appropriate data to the structure 
    trc_struct.NumFrames = size(M,1);
    trc_struct.DataRate = rate;
    
    trc_struct.data.(fields{1}) = M(:,1);
    trc_struct.data.(fields{2}) = M(:,2);
    
    for col=3:numel(fields);
       marker_data = M(:,(col-2)*3:(col-2)*3+2);
       trc_struct.data.(fields{col}) = marker_data;
    end
end