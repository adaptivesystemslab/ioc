function returnVal = lookupTableDumas(parameter, segment, gender, lengthSegment, massSegment)
% https://www.sciencedirect.com/science/article/pii/S0021929006000728#tbl2

% length [%]: normalized to height of the participant. need to multiply by the
%   participant height
% mass [%]: normalized weight of the participant. need to multiply by the
%   participant weight
% com [m]
% inertial []

dumasTable = readtable('Dumas2006.csv');  

findSegment = find(strcmpi(segment, table2array(dumasTable(:, 1))));
findGender = find(strcmpi(gender, table2array(dumasTable(:, 4))));
row = intersect(findSegment, findGender);

switch parameter
    case 'length'
        ind = 5;
        switch gender
            case 'm'
                avgHeight = 1.77;
                
            case 'f'
                avgHeight = 1.61;
        end
        
        returnVal = table2array(dumasTable(row, ind)) / 1000; % [mm] to [m]
        returnVal = returnVal / avgHeight; 
        
    case 'mass'
        ind = 6;
        returnVal = table2array(dumasTable(row, ind)) / 100; % [%]
        
    case 'com'
        ind = 7:9;
        outputRaw = table2array(dumasTable(row, ind))' / 100; % [%]
        
        returnVal = outputRaw * lengthSegment;
         
    case 'inertial'
        ind = 10:15;
        outputRaw = table2array(dumasTable(row, ind))';
        outputRaw2 = [...
            outputRaw(1) outputRaw(4) outputRaw(5)
            outputRaw(4) outputRaw(2) outputRaw(6)
            outputRaw(5) outputRaw(6) outputRaw(3)] / 100;
        
        returnVal = ((outputRaw2 * lengthSegment) .^ 2) * massSegment;
end
end