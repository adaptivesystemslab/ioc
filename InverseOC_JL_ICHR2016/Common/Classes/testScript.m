% test class structure
function testScript
    % preamble
%     basePath = 'D:\MyDocument\MotionData\Lowerbody_healthy1_2011-11\';
    basePath = 'D:\MyDocument\MotionData\Lowerbody_healthy2_2013-07\';

    testExerciseStruct([basePath 'Subject5\Session1\KFEO_STD\'], 'Healthy2');
%     testIndividualDH([basePath 'Subject5\Session1\KFEO_STD\']);
end

function testExerciseStruct(basePath, project)
    exerciseHandle = exerciseDataHandle(basePath);
end

function testIndividualDH(basePath)
    % loading Cortex data
%     headerPath = [basePath 'Cortex\MocapData_Cortex.header'];
%     dataPath = [basePath 'Cortex\MocapData_Cortex.data'];
%     mocapData = arsLoader(headerPath, dataPath);
    
    % loading IMU data
    headerPath = [basePath 'Shimmer\SensorData_Waist.header'];
    dataPath = [basePath 'Shimmer\SensorData_Waist.data'];
    imuData{1} = arsLoader(headerPath, dataPath);
    
    headerPath = [basePath 'Shimmer\SensorData_RKnee.header'];
    dataPath = [basePath 'Shimmer\SensorData_RKnee.data'];
    imuData{2} = arsLoader(headerPath, dataPath);
    
    headerPath = [basePath 'Shimmer\SensorData_RAnkle.header'];
    dataPath = [basePath 'Shimmer\SensorData_RAnkle.data'];
    imuData{3} = arsLoader(headerPath, dataPath);
    
    % loading EKF data
    headerPath = [basePath 'EKF\2014_02_13\ekf.header'];
    dataPath = [basePath 'EKF\2014_02_13\ekf.data'];
    ekfData = arsLoader(headerPath, dataPath);
%     
%     % loading manual segment data
    headerPath = [basePath 'Segmentation_manual\SegmentData_Manual_Manual.header'];
    dataPath = [basePath 'Segmentation_manual\SegmentData_Manual_Manual.data'];
    segManData = arsLoader(headerPath, dataPath);
%     
    clc
    shimmerData = imuData{1}
    ekfData
    segManData
end