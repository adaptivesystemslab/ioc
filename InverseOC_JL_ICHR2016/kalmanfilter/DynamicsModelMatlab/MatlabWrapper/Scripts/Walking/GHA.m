function [R, gyro_bias] = GHA(imu_static,imu_dynamic,grav_axis,m_axis,...
    init_rot_axis,static_indxs,dynamic_indxs)

%Defaults 
if ~exist('grav_axis','var')
   grav_axis = 3; 
end
if ~exist('m_axis','var')
   m_axis = 2; 
end
if ~exist('static_indxs','var')
   static_indxs = 1000:(size(imu_static.accelerometerCalibrated,1)-1000); 
end
if ~exist('dynamic_indxs','var')
   dynamic_indxs = 1000:(size(imu_dynamic.gyroscopeCalibrated,1)-1000);
end

%Static Pose Learning Rate
nu = 0.05;
R_accel = zeros(3,3);

accel_still = imu_static.accelerometerCalibrated(static_indxs,:);
gyro_bias = mean(imu_static.gyroscopeCalibrated(static_indxs,:));

Ax = zeros(3,1);
Ax(grav_axis) = 1;

w_axis = m_axis;
Px = zeros(3,1);
Px(m_axis) = 1;

Lx = zeros(3,1);
last_axis = find(Ax+Px == 0,1);
Lx(last_axis) = 1;

for i = 1:size(accel_still,1)
    x = accel_still(i,:)';
    x = x/norm(x);

    deltaAx = dot(Ax,x)*(x);

    Ax = Ax + nu*deltaAx;
    Ax = Ax/norm(Ax);

    %Do the cross products
    %Y is in grav Z is perp
    if(w_axis == 3 && last_axis == 1)
        Lx = cross(Ax,Px);
        Px = cross(Lx,Ax);
    else
        %Z is up X is zero G
        Lx = cross(Ax,Px);
        Lx = Lx/norm(Lx);
        Px = cross(Lx,Ax);
        Px = Px/norm(Px);
    end
    R = zeros(3,3);
    R(:,grav_axis) = Ax;
    R(:,m_axis) = Px;
    R(:,last_axis) = Lx;
end

% mean_accel = mean((R'*accel_still')');
% if mean_accel(grav_axis) < 0
%     R = R*rotx(pi);
% end

R_accel(:,:) = R;

%Learning Rate
nu = 0.02;

%Final rotatino matrices
R = R_accel;

mes_dynamic = imu_dynamic.gyroscopeCalibrated(dynamic_indxs,:);
mes_dynamic = mes_dynamic - repmat(gyro_bias,size(mes_dynamic,1),1);

%Get axes based on accel calibrated data
R = R_accel;
Ax = R(:,grav_axis);
%Px = R(:,w_axis);
Px = init_rot_axis;
Lx = R(:,last_axis);

for i = 1:size(mes_dynamic,1)
    x = mes_dynamic(i,:)';
    %Subtract off motion in gravity axis
    x = x - dot(x,Ax)*Ax;
    deltaPx = dot(Px,x)*(x);
    deltaLx = dot(Lx,x)*(x - dot(x,Lx)*Lx);
    deltaAx = dot(Ax,x)*(x -dot(x,Px)*Px -dot(x,Lx)*Lx);
    
    Px = Px + nu*deltaPx;
    Px = Px/norm(Px);
    
    
    Lx = cross(Ax,Px);
    Lx = Lx/norm(Lx);
    Px = cross(Lx,Ax);
    Px = Px/norm(Px);
    
    R = zeros(3,3);
    R(:,grav_axis) = Ax;
    R(:,m_axis) = Px;
    R(:,last_axis) = Lx;
end

%Resulting rotation to rotate IMU data to desired frame definition
R  = R';

end

