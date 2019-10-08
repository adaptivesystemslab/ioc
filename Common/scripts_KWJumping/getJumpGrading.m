function [targetLengths,jumpGrades] = getJumpGrading(partNum)
% gender = 'M' or 'F'
% height = integer, in [cm]
% weight = integer, in [lbs]
% S2sc = shoulder marker to shoulder center vertical distance, in [cm]

% Read in "P#_jump_grading.csv" file
mydir = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)); 
csvPath = [newdir 'data\P' partNum '_filtered/P' partNum '_jump_grading.csv'];
if(exist(csvPath,'file')~=2)
    targetLengths = 0;
    jumpGrades = 'empty';
else
    fid = fopen(csvPath);
    C = textscan(fid, '%d%s%s%s%s%s%s', 'Delimiter', ',', 'HeaderLines', 2);
    fclose(fid);

    targetLengths = [C{1,1}];
    % jumpGrades ordered as on data collection notes: [target_1,set_1],[T2,S1],[T3,S1],[T1,S2],[T2,S2][T3,S2]
    for i = 1:numel(C{1,1}) % num targets jumped to
        for j = 1:6 % num jumps at each target
            jumpGrades{i,j} = C{1,j+1}(i);
        end
    end
    
    % put jumpGrades in prefered order: [target_1,set_1],[T1,S2],[T2,S1],[T2,S2],[T3,S1][T3,S2]
    jumpGradesTmp = jumpGrades;
    jumpGradesTmp(2,:) = jumpGrades(4,:);
    jumpGradesTmp(3,:) = jumpGrades(2,:);
    jumpGradesTmp(4,:) = jumpGrades(5,:);
    jumpGradesTmp(5,:) = jumpGrades(3,:);
    
    % vectorize jumpGrades into "jump ordering", as described in Scientific
    % Data paper submission
    jumpGrades = reshape(jumpGradesTmp',36,1);
    
end
