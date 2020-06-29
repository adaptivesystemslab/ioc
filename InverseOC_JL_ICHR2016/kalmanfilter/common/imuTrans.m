function transform = imuTrans(name, frame, T_f_raw_svdAvg, T_f_0_svdAvg, T_f_imu_svdAvg)
    %Transformation From IMU frame into World Frame
    switch nargin
        case 2
            transform.name = name;
            transform.frame = frame.name;
            transform.T_f_imu = [];
            transform.T_f_0 = [];
            transform.R = [];
            transform.t = [];
            
        otherwise
            %Add Sensor to Model
            transform.name = name;
            transform.frame = frame.name;
            transform.T_f_imu = T_f_raw_svdAvg;
            transform.T_f_0 = T_f_0_svdAvg;
            transform.R = T_f_imu_svdAvg(1:3, 1:3);
            transform.t = T_f_imu_svdAvg(1:3, 4);
    end
end