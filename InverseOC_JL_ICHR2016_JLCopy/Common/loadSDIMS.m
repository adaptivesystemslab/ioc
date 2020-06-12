function model = loadSDIMS(file_sdims)
%Go through the file line by line creating the model 

%This is a tree structure which will keep track of the model, later it is
%converted to XML doc readable by our EKF. The root node only contains the
%world frame information.

%Create empty tree
model = tree();
%Keeps track of current node
current_node = 0;

fid = fopen(file_sdims);

tline = fgets(fid);
while ischar(tline)
    %Skip empty lines
    if strcmp(sprintf('%x',tline),'da')
        tline = fgets(fid);
        continue;
    end
    
    %Model is empty so we just started reading the file and we found world acceleration data
    if isempty(model.Node{1}) && ~isempty(strfind(tline,'LinAcceleration'))
        g = regexp(tline, '\d*\.?\d*','Match');
        g = str2double(g);
        node.name = 'world';
        node.g = g;
        model = model.set(1,node);
        current_node = 1;
    end
    
    %We found a new body to add and already dealt with world
    if ~isempty(strfind(tline,'Link')) && current_node ~= 0
        %Get the name of the link
        tline = regexprep(tline,'\r\n|\n|\r','');
        strings = strsplit(tline,' ');
        node_name = strings{2};
        
        node = struct();
        node.name = node_name;
        %Really this should be recursive implementation but i'm lazy
        %Basically read until we hit empty line decplaring end of node
        
        %Get the data
        tline = fgets(fid);
        while ~isempty(regexprep(tline,'\r\n|\n|\r',''))
            
            %Put into struct 
            strings = strsplit(tline,' ');
            node.(strtrim(strings{1})) = strtrim(strings(2:end));
            tline = fgets(fid);
        end
        %Add the node to the model and mode on
        [model, current_node] = model.addnode(current_node,node);
    end
    
    %Get the next line of the file
    tline = fgets(fid);
    
    %Check if it is the end of a branch
    if strfind(tline,'}')
        %It is the end of the branch so we go back to the parents
       current_node= model.getparent(current_node);
    end
    
end
