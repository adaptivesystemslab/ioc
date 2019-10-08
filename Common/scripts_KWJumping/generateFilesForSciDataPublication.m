% Used to create Participant Data .CSV Files for Nature Scientific Data
% Publication

% EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\CloudSLAM\lg_persontracking\';
% addpath('..\..\');
EKFCodePath = 'C:\Users\kgwester\Documents\ResearchWork\MatlabWrapper\';
addpath(EKFCodePath);


savePath = ['C:\Users\kgwester\Documents\ResearchWork\Jumping Dataset Publication - Nature Scientific Data\Data\'];


% Choose files to generate and save
jump_grading = 1;
model_links = 1;
world2base = 1;
takeoff_land = 1;
shift_align_values = 1;
joint_angle_data = 0;
Rshldr_Rhip_transforms = 0;


partsToRun = {'02','03','04','05','06','07','08','09','10','11','12',...
              '13','14','15','16','17','18','19','20','21','22'};
% partsToRun = {'04'};


mydir = pwd;
idcs = strfind(mydir,'\');
newdir = mydir(1:idcs(end));


for i_part = 1:numel(partsToRun)
    
    partNum = partsToRun{i_part};
%     [age,gender,height,weight,S2sc,calibPose] = getPartData(partNum);
%     S2sc = S2sc/100;

    dataFolder = ['P' partNum '_filtered/'];
    targetDist = {'55','70','85'};


    loadFilePath = [newdir 'results\JA_P' partNum '.mat'];
    if(exist(loadFilePath,'file')~=2)
        disp(['JA_P' partNum '.mat NOT FOUND']);
    else
        load(loadFilePath);
        
        savePath_part = [savePath 'P' partNum '_Data\'];
        
        
        if(jump_grading)
            fileName = ['P' partNum '_jump_grades.csv'];
%             jump_grades = reshape(JA.jumpGrades',36,1);
            jump_grades = JA.jumpGrades; % now in "jump ordering" in JA struct, [36 x 1] vector
            jumpNum = [1:36]';
            header = 'Jump Number, Jump Grade';
            
            fid = fopen([savePath_part fileName],'w');
            fprintf(fid,[header '\n']);
            for line = 1:36
                fprintf(fid,'%d,%s\n',jumpNum(line),char(jump_grades{line}));
            end
            fclose(fid);
            % Can't use "saveCSVWithHeader" b/c of non-numeric characters
        end
        

        if(model_links)
            fileName = ['P' partNum '_model_links.csv'];
            modelLinkNames = fieldnames(JA.modelLinks);
            header = 'Model Link Name, Global X [m], Global Y [m], Global Z [m]';
            
            fid = fopen([savePath_part fileName],'w');
            fprintf(fid,[header '\n']);
            for line = 1:numel(modelLinkNames)
                if(numel(JA.modelLinks.(modelLinkNames{line})) == 1)
                    fprintf(fid,'%s,%f,%f,%f\n',modelLinkNames{line},0,0,JA.modelLinks.(modelLinkNames{line}));
                else % 3-element vector
                    fprintf(fid,'%s,%f,%f,%f,\n',modelLinkNames{line},JA.modelLinks.(modelLinkNames{line}));
                end
            end
            fclose(fid);
        end
        
        
        if(world2base)
            fileName = ['P' partNum '_world2base.csv'];
            data = [[1:36]', reshape(JA.world2base,36,16)];
            header = 'Jump Number, Vectorized World2Base Transform';
            saveCSVWithHeader([savePath_part fileName],header,data);
        end
        
        
        if(takeoff_land)
            fileName = ['P' partNum '_takeoff_land_frames.csv'];
            data = [[1:36]', reshape(JA.TOFrame,36,1), reshape(JA.LandFrame,36,1)];
            header = 'Jump Number, Takeoff Frame, Land Frame';
            saveCSVWithHeader([savePath_part fileName],header,data);
        end
        
        
        if(shift_align_values)
            land_frames_align = [JA.LandFrame(:,1)-JA.LandFrame(12,1),...
                JA.LandFrame(:,2)-JA.LandFrame(12,2),...
                JA.LandFrame(:,3)-JA.LandFrame(12,3)];
            shift_frames_align = JA.LandFrame_shiftValues;
            fileName = ['P' partNum '_shift_align_values.csv'];
            data = [[1:36]', reshape( (land_frames_align + shift_frames_align),36,1)];
            % NOTE: for land_frames_align and shift_frames_align, positive 
            % value means shift EARLIER in time (remove frames from front 
            % of trajectory), while negative value means shift LATER in 
            % time (add placeholder frames to front of trajectory)
            header = 'Jump Number, Frame Shift Value';
            saveCSVWithHeader([savePath_part fileName],header,data);
        end
        
        
        
        % Generate and save files for individual jump recordings
        for i_targ = 1:3
            for i_set = 1:2
                for i_jump = 1:6

                    currJump = 12*(i_targ-1) + 6*(i_set-1) + i_jump;
                    
                    
                    if(joint_angle_data)
                        JA_data = JA.targ(i_targ).jump(6*(i_set-1) + i_jump).data;
                        fileName = ['P' partNum '_JA_jump_' num2str(currJump) '.csv'];
                        
                        header = [JA.jointNames(1:end-1); repmat({','},1,size(JA.jointNames,2)-1)]; %insert commas
                        header = header(:)';
                        header = cell2mat(header);
                        header = strcat('Frame,',header,JA.jointNames{end}); % so there is no delimeter at the end of the header line
                        
                        data = [[1:size(JA_data,1)]', JA_data];
                        saveCSVWithHeader([savePath_part 'Joint Angles\' fileName],header,data);
                    end
                    
                    
                    if(Rshldr_Rhip_transforms)
                        R_matNames = {'shldr_L','shldr_R','hip_L','hip_R'};
                        
                        for i_R = 1:4
                            R_mat = JA.targ_R_mat(i_targ).jump(6*(i_set-1) + i_jump).(R_matNames{i_R});
                            data = [[1:size(R_mat,1)]',reshape(R_mat,size(R_mat,1),9)];
                            fileName = ['P' partNum '_JA_jump_' num2str(currJump) '_R-' R_matNames{i_R} '.csv'];
                            header = ['Frame, Vectorized R-' R_matNames{i_R} ' matrix'];
                            saveCSVWithHeader([savePath_part 'Joint Angles\' fileName],header,data);
                        end
                    end
                    
                    
                    
                    
                end
            end
        end
        
        
        
        disp(['Finished generating files for participant ' partNum]);
        
    end
end

                    
                    
                    
                    
                    
                    
                    