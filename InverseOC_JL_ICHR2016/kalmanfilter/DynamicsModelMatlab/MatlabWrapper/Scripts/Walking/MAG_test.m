%% Calibrate Mag

addpath(genpath('C:\asl_git\kalmanfilter\ik_framework\common\ars'));
addpath('../../');

imu_path = 'C:\aslab\data\Gait_lowerbody_2019\Pilot Data\magtest\imu\';
imu = arsLoader([imu_path '0006667D7185_2019_06_25_11_40_37.header']);

mag = [imu.data.MagnetometerXUncalibrated imu.data.MagnetometerYUncalibrated ...
    imu.data.MagnetometerZUncalibrated];

xy = [mag(:,1) mag(:,2)];
xz = [mag(:,1) mag(:,3)];
yz = [mag(:,2) mag(:,3)];

%Get rid of outliers by looking at distance to mean
xy_mean = mean(xy);
xy_mean_dist = sqrt(sum((xy-repmat(xy_mean,size(xy,1),1)).^2,2));
xy_fixed = xy(xy_mean_dist < 150,:);
xy_xmid = min(xy_fixed(:,1))+(max(xy_fixed(:,1))-min(xy_fixed(:,1)))/2;
xy_ymid = min(xy_fixed(:,2))+(max(xy_fixed(:,2))-min(xy_fixed(:,2)))/2;
xy_fixed = xy_fixed -repmat([xy_xmid xy_ymid],size(xy_fixed,1),1);

xz_mean = mean(xz);
xz_mean_dist = sqrt(sum((xz-repmat(xz_mean,size(xz,1),1)).^2,2));
xz_fixed = xz(xz_mean_dist < 150,:);
xz_xmid = min(xz_fixed(:,1))+(max(xz_fixed(:,1))-min(xz_fixed(:,1)))/2;
xz_zmid = min(xz_fixed(:,2))+(max(xz_fixed(:,2))-min(xz_fixed(:,2)))/2;
xz_fixed = xz_fixed -repmat([xz_xmid xz_zmid],size(xz_fixed,1),1);

yz_mean = mean(yz);
yz_mean_dist = sqrt(sum((yz-repmat(yz_mean,size(yz,1),1)).^2,2));
yz_fixed = yz(yz_mean_dist < 150,:);
yz_ymid = min(yz_fixed(:,1))+(max(yz_fixed(:,1))-min(yz_fixed(:,1)))/2;
yz_zmid = min(yz_fixed(:,2))+(max(yz_fixed(:,2))-min(yz_fixed(:,2)))/2;
yz_fixed = yz_fixed -repmat([yz_ymid yz_zmid],size(yz_fixed,1),1);

%Now PCA the results and scale 
x_offset = (xy_xmid + xz_xmid)/2;
y_offset = (xy_ymid + yz_ymid)/2;
z_offset = (xz_zmid + yz_zmid)/2;

mag_bias = [x_offset y_offset z_offset];
mag_ub = mag - repmat(mag_bias,size(mag,1),1);

mag_scalex = (max(xz_fixed(:,1)) - min(xz_fixed(:,1)))/2;
mag_scaley = (max(yz_fixed(:,1)) - min(yz_fixed(:,1)))/2;
mag_scalez = (max(xz_fixed(:,2)) - min(xz_fixed(:,2)))/2;
mag_scale = [mag_scalex mag_scaley mag_scalez];
mag_ub_scale = mag_ub./repmat(mag_scale,size(mag_ub,1),1);


figure(1);clf;hold on;
plot(xy_fixed(:,1),xy_fixed(:,2),'r*');
plot(xz_fixed(:,1),xz_fixed(:,2),'g*');
plot(yz_fixed(:,1),yz_fixed(:,2),'b*');


figure(2);clf;hold on;
plot3(mag_ub(:,1),mag_ub(:,2),mag_ub(:,3),'.');
grid on;

figure(2);clf;hold on;
plot3(mag_ub_scale(:,1),mag_ub_scale(:,2),mag_ub_scale(:,3),'.');
grid on;

%% Try to apply this to the outside walking IMU 

%imu_path = 'C:\aslab\data\Gait_lowerbody_2019\Data Collection\Subject04_2019_06_07\imu\';
imu_path = 'C:\aslab\data\Gait_lowerbody_2019\MagCalib\';
imu = arsLoader([imu_path '0006667D7185_2019_06_17_14_39_47.header']);

mag = [imu.data.MagnetometerXUncalibrated imu.data.MagnetometerYUncalibrated ...
    imu.data.MagnetometerZUncalibrated];
mag_ub = mag - repmat([x_offset y_offset z_offset],size(mag,1),1);


mag_ub_scale = mag_ub./repmat(sqrt(sum(mag_ub.^2,2)),1,3);


r = vrrotvec(mean(mag_ub_scale(1:100,:)),[0 0 1]);
R = vrrotvec2mat(r);
%mag_ub_rot_scale = (R*mag_ub_rot_scale')'

figure(3);clf;hold on;
plot(mag_ub_scale(:,1),'r');
plot(mag_ub_scale(:,2),'g');
plot(mag_ub_scale(:,3),'b');

% Project onto XY and figure out angle
mag_xy = mag_ub_scale(:,[1 2]);
mag_xy_scale = mag_xy./repmat(sqrt(sum(mag_xy.^2,2)),1,2);
mag_xy_init = mean(mag_xy_scale(1:100,:));
angles = zeros(size(mag_xy,1),1); 
for i=1:size(mag_xy,1)
    angles(i) = acos(dot(mag_xy_init,mag_xy_scale(i,:)));
end

figure(4)
plot(rad2deg(angles));

%%

figure(4);clf;hold on;axis([-2 2 -2 2 -2 2]);grid on
hn=arrow3([0 0 0],[mag_ub_scale(1,:)]);
for i=1300:size(mag_ub_scale)
    delete(hn)
    hn=arrow3([0 0 0],[mag_ub_scale(i,:)]);
    drawnow;
    pause(0.01);
end