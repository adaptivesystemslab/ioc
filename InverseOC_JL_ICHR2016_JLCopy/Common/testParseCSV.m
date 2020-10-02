% Testing parseCSV function
basePath1 = 'C:\Documents\aslab\data\Lowerbody_healthy1_2011-11\Subject30\Session1\SQUA_STD_SLO1\Cortex\';
basePath2 = 'C:\Documents\aslab\data\Lowerbody_healthy1_2011-11\Subject31\Session1\SQUA_STD_SLO1\Cortex\';
basePath3 = 'C:\Documents\aslab\data\Lowerbody_healthy1_2011-11\Subject21\Session1\SQUA_STD_SLO1\Cortex\';
head = 'MocapData_Cortex.header';
dat = 'MocapData_Cortex.data';

basePath4 = 'C:\Users\andre\Documents\MATLAB\Doppel';

%% Parse into Matlab
cd(basePath1);
squats1to3 = parseCSV(dat);

cd(basePath2);
squats4to6 = parseCSV(dat);

cd(basePath4);
euler1to3 = parseCSV('1-3Euler_MocapData_Cortex.data');

all = parseCSV('all_MocapData_Cortex.data');

dataset =  parseCSV('reformatted.csv');

% cd(basePath3);
% all = parseCSV(data);

%% Try plotting the data that was read in
plotHip = 0;
plotKnee = 1;
data = squats1to3;
data1 = squats4to6;
euler = euler1to3;
%%
figure (6), hold on
plot(all.KneeRight_x)
plot(all.KneeRight_y)
plot(all.KneeRight_z)
legend('x','y','z')
%%

if(plotHip)
    figure (10)
    subplot(3,1,1), hold on;
    plot(data.HipRight_qx);
    plot(data.HipRight_qy);
    plot(data.HipRight_qz);
    plot(data.HipRight_qw);
    title('HIPS')        
    legend('qx', 'qy', 'qz', 'qw')
    xlim([0 150])
    
    subplot(3,1,2), hold on
    plot(euler.HipRightAngle_x);
    plot(euler.HipRightAngle_y);
    plot(euler.HipRightAngle_z);
    legend('x','y','z')
    xlim([0 150])
    
    subplot(3,1,3), hold on
    plot(all.HipRight_x(135:285,:))
    plot(all.HipRight_y(135:285,:))
    legend('x','y')
    xlim([0 150])
    
    
end

if(plotKnee)
    subplot(2,1,1), hold on;
    plot(data.KneeRight_qx);
    plot(data.KneeRight_qy);
    plot(data.KneeRight_qz);
    plot(data.KneeRight_qw);
    title('QUATERNIONS R KNEE')        
    legend('qx', 'qy', 'qz', 'qw')
    xlim([0 150])
   
%     % Recalculate Euler angles
%     n = length(data.KneeRight_qx);
%     z = zeros(1,n);
%     y = zeros(1,n);
%     x = zeros(1,n);
%     for i = 1:n
%         quat = [data.KneeRight_qx(i), data.KneeRight_qy(i), data.KneeRight_qz(i), data.KneeRight_qw(i)];
%         eul = quat2euler(quat);
%         x(i) = eul(1);
%         y(i) = eul(2);
%         z(i) = eul(3);
%     end
%     
%     subplot(3,1,2), hold on
%     plot(euler.KneeRightAngle_x);
%     plot(euler.KneeRightAngle_y);
%     plot(euler.KneeRightAngle_z);
    
%     plot(x)
%     plot(y)
%     plot(z)
%     
%     title('EULER ANGLES R KNEE')     
%     legend('x','y','z')
%     xlim([0 150])
%     ylabel('Degrees');
    
    subplot(2,1,2), hold on
    plot(all.HipRight_y(135:285,:))
    legend('y')
    xlim([0 150])
    title('Y POSITION R HIP');

end

Hip = 1;
if(Hip)
    figure(9)
    subplot(2,1,1), hold on;
    plot(data.HipRight_qx);
    plot(data.HipRight_qy);
    plot(data.HipRight_qz);
    plot(data.HipRight_qw);
    title('QUATERNIONS R HIP')        
    legend('qx', 'qy', 'qz', 'qw')
    xlim([0 150])
   
%     % Recalculate Euler angles
%     n = length(data.HipRight_qx);
%     z = zeros(1,n);
%     y = zeros(1,n);
%     x = zeros(1,n);
%     for i = 1:n
%         quat = [data.HipRight_qx(i), data.HipRight_qy(i), data.HipRight_qz(i), data.HipRight_qw(i)];
%         eul = quat2euler(quat);
%         x(i) = eul(1);
%         y(i) = eul(2);
%         z(i) = eul(3);
%     end
%     
%     subplot(3,1,2), hold on
% %     plot(euler.KneeRightAngle_x);
% %     plot(euler.KneeRightAngle_y);
% %     plot(euler.KneeRightAngle_z);
%     
%     plot(x)
%     plot(y)
%    plot(z)
%     
%     title('EULER ANGLES R HIP')     
%     legend('x','y','z')
%     xlim([0 150])
%     ylabel('Degrees');
    
    subplot(2,1,2), hold on
    plot(all.HipRight_y(135:285,:))
    legend('y')
    xlim([0 150])
    title('Y POSITION R HIP');

end


ankle = 1;

if(ankle)
    figure(8)
    subplot(2,1,1), hold on;
    plot(data.AnkleRight_qx);
    plot(data.AnkleRight_qy);
    plot(data.AnkleRight_qz);
    plot(data.AnkleRight_qw);
    title('QUATERNIONS R ANKLE')        
    legend('qx', 'qy', 'qz', 'qw')
    xlim([0 150])
%    
%     % Recalculate Euler angles
%     n = length(data.AnkleRight_qx);
%     z = zeros(1,n);
%     y = zeros(1,n);
%     x = zeros(1,n);
%     for i = 1:n
%         quat = [data.AnkleRight_qx(i), data.AnkleRight_qy(i), data.AnkleRight_qz(i), data.AnkleRight_qw(i)];
%         eul = quat2euler(quat);
%         x(i) = eul(1);
%         y(i) = eul(2);
%         z(i) = eul(3);
%     end
%     
%     subplot(3,1,2), hold on
% %     plot(euler.KneeRightAngle_x);
% %     plot(euler.KneeRightAngle_y);
% %     plot(euler.KneeRightAngle_z);
%     
%     plot(x)
%     plot(y)
%     plot(z)
%     
%     title('EULER ANGLES R ANKLE')     
%     legend('x','y','z')
%     xlim([0 150])
%     ylabel('Degrees');
    
    subplot(2,1,2), hold on
    plot(all.HipRight_y(135:285,:))
    legend('y')
    xlim([0 150])
    title('Y POSITION R HIP');

end
