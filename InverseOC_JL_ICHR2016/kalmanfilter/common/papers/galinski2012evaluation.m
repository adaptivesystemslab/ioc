function transform = galinski2012evaluation(dataInstance, modelInstance)

end
% 
% clearvars
% addpath(genpath('..\2018_06_08')); % add 'common' to the filepath
% if currFileEntry.fileId == 'HAAO_STD_NON1' 
%     filepath = currFileEntry.filePathImuKneeRight;
%     imudata = arsLoader(filepath);
%     gyro = imudata.gyroscopeCalibrated;
%     
%     for i=1:length(gyro)
%         for j=1:3
%             norm_x(i)=sqrt(gyro(i,j)^2+gyro(i,j)^2+gyro(i,j)^2);
%         end
%     end
% 
% [Max_X,I] = max(norm_x);
% X_cal=[gyro(I,1),gyro(I,2),gyro(I,3)]/(sqrt(gyro(I,1)^2+gyro(I,2)^2+gyro(I,3)^2));
% end
% 
% if currFileEntry.fileId == 'HEFO_STD_NON1' 
%     filepath = currFileEntry.filePathImuKneeRight;
%     imudata = arsLoader(filepath);
%     gyro = imudata.gyroscopeCalibrated;
%     
%     for s=1:length(gyro)
%         for t=1:3
%             norm_y(s)=sqrt(gyro(s,t)^2+gyro(s,t)^2+gyro(s,t)^2);
%         end
%     end
% 
% [Max_Z,R] = max(norm_y);
% Z_cal=[gyro(R,1),gyro(R,2),gyro(R,3)]/(sqrt(gyro(R,1)^2+gyro(R,2)^2+gyro(R,3)^2));
% end
% 
% Y_cal = cross(X_cal,Z_cal);
% 
% Rot_Struct.name = imudata.name;
% Rot_Struct.R_calMatrix = [X_cal',Y_cal',Z_cal']; 
% Rot_Struct.t_calTrans = [];
% 
% disp(Rot_Struct.R_calMatrix);
