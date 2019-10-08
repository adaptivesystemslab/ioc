% load JA data from multiple participants, get all jump grades from all
% jumps into single array



partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};

jumpGradeTypes = {'B','SB','P','P*','SF','F'};
jumpGradesAll = {};
jumpGradesAll_num = [];

for partNum = 1:numel(partsToRun)
    % load JA data struct from certain participant
    mydir = pwd;
    idcs = strfind(mydir,'\');
    newdir = mydir(1:idcs(end)-1); 
    loadFilePath = [newdir '\results\JA_P' partsToRun{partNum} '.mat'];
    if(exist(loadFilePath,'file')~=2)
        disp(['JA_P' partNum '.mat NOT FOUND']);
    else
        load(loadFilePath); % JA data struct, human data
        
        for jumpNum = 1:36
            jumpGradesAll{jumpNum,partNum} = JA.jumpGrades{jumpNum};
            
%             gradeTypeNum = find(strcmp(JA.jumpGrades{jumpNum},jumpGradeTypes));
            jumpGradesAll_num(jumpNum,partNum) = find(strcmp(JA.jumpGrades{jumpNum},jumpGradeTypes));
            
            % Change grade numbering from 1 to 5, with P* = 3.3 (for excel color scale formatting)
%             if(jumpGradesAll_num(jumpNum,partNum) == 4)
%                 jumpGradesAll_num(jumpNum,partNum) = 3.3;
%             elseif(jumpGradesAll_num(jumpNum,partNum) >4)
%                 jumpGradesAll_num(jumpNum,partNum) = jumpGradesAll_num(jumpNum,partNum) - 1;
%             end
            
        end
        
    end
end

