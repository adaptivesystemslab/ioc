function returnVal = lookupTableDumas(parameter, segment, gender, lengthSegment, massSegment)
% https://www.sciencedirect.com/science/article/pii/S0021929006000728#tbl2

% segment = 'forearm';
% gender = 'm';
% parameter = 'inertial';

dumasTable = readtable('Dumas2006.csv');  

findSegment = find(strcmpi(segment, table2array(dumasTable(:, 1))));
findGender = find(strcmpi(gender, table2array(dumasTable(:, 4))));
row = intersect(findSegment, findGender);

switch parameter
    case 'mass'
        ind = 6;
        returnVal = table2array(dumasTable(row, ind)) / 100;
        
    case 'com'
        ind = 7:9;
        outputRaw = table2array(dumasTable(row, ind))' / 100;
        
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

% switch parameter
%     case 'mass'
%         switch segment
%             case 'foot' % AJCtoMH15
%                 returnVal = 1.2/100;
%                 
%             case 'leg' % KJCtoAJC
%                 returnVal = 4.8/100;
%                 
%             case 'thigh' % HJCtoKJC
%                 returnVal = 12.3/100;
%                 
%             case 'pelvis' % LJCtoHJC
%                 returnVal = 14.2/100;
%                 
%             case 'torso' % CJCtoLJC
%                 returnVal = 33.3/100;
%                 
%             case 'arm' % SJCtoEJC
%                 returnVal = 2.4/100;
%                 
%             case 'forearm' % EJCtoWJC
%                 returnVal = 1.7/100;
%                 
%             case 'hand' % WJCtoMH25
%                 returnVal = 0.6/100;
%                 
%             case 'headneck' % CJCtoHV
%                 returnVal = 6.7/100;
%         end
%         
%     case 'com'
%         
%     case 'inertial'
% end