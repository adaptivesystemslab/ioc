function [age,gender,height,weight,S2sc,calibPose] = getPartData(partNum)
% gender = 'M' or 'F'
% height = integer, in [cm]
% weight = integer, in [lbs]
% S2sc = shoulder marker to shoulder center vertical distance, in [cm]

if(ischar(partNum) || isstring(partNum))
    partNum = str2num(partNum);
end

% Read in "additional_participant_data.csv" file
mydir = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)); 
csvPath = [newdir 'data\additional_participant_data.csv'];
fid = fopen(csvPath);
C = textscan(fid, '%d%d%c%d%d%f%c%d', 'Delimiter', ',', 'HeaderLines', 1);
fclose(fid);

age = C{1,2}(partNum);
gender = C{1,3}(partNum);
height = C{1,4}(partNum);
weight = C{1,5}(partNum);
S2sc = C{1,6}(partNum); 
calibPose = C{1,7}(partNum);
% physicalActivity = C{1,8}(partNum); % participant;s estimate of how many hours per week they are physically active
