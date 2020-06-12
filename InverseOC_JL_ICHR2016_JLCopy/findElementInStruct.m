function output = findElementInStruct(baseStruct, query)

    output = struct;
    
    for i=1:length(baseStruct)
       if strcmp(baseStruct(i).name, query)
           output = baseStruct(i);
           return;
       end
    end
    
end

