function subjectInfo = loadPatientData(fileStack, shimmerDevicesFiles, ekfSubpath)
    % load all patient data. ekf, seg, if exists. 
    shimmerPathPrefix = 'Shimmer';
    cortexPathPrefix = 'Cortex';
    
    % path prefix define
    ekfPathPrefix = 'EKF';
    
    subjectInfoInd = 0;
    
    subjectInfo = [];
    
    if ~iscell(fileStack)
        % the function expects a cell
        fileStackTemp = fileStack;
        clear fileStack;
        fileStack{1}.subjectNumber = fileStackTemp.subjectNumber;
        fileStack{1}.sessionNumber = fileStackTemp.sessionNumber;
        fileStack{1}.exerciseName = fileStackTemp.exerciseName;
        fileStack{1}.exerciseType = fileStackTemp.exerciseType;
        fileStack{1}.filePath = fileStackTemp.filePath;
    end
    
    for ind_fileStack = 1:length(fileStack)
        currSubjNum = fileStack{ind_fileStack}.subjectNumber;
        currSessNum = fileStack{ind_fileStack}.sessionNumber;
        currExer = fileStack{ind_fileStack}.exerciseName;
        currExerPath = fileStack{ind_fileStack}.filePath;

        fprintf('Currently loading: %s\n', currExerPath);
        
        % defining filepaths
        shimmerPath = fullfile(currExerPath, shimmerPathPrefix);
        ekfPathTemp = fullfile(currExerPath, ekfPathPrefix);
        cortexPath = fullfile(currExerPath, cortexPathPrefix, 'MocapData_Cortex.data');
        
        if ~exist('ekfSubpath', 'var') || isempty(ekfSubpath)
            % need to determine the EKF subfolder name (assume the last
            % entry is the one we want), if no specific subfolder is passed
            % in
            ekfSubfolderDir = dir(ekfPathTemp);
            
            if isempty(ekfSubfolderDir)
                continue
            else
                ekfPath = fullfile(ekfPathTemp, ekfSubfolderDir(end).name);
            end
        else
            ekfPath = fullfile(ekfPathTemp, ekfSubpath);
        end
        
        % now need to parse the contents of the EKF folder
        % pull out joint angle data
        ekfFilefindingDir = dir(fullfile(ekfPath, 'JointAngle*.data'));
        if ~isempty(ekfFilefindingDir)
            % separated layout
            jointDataPath = fullfile(ekfPath, ekfFilefindingDir.name);
        else
            ekfFilefindingDir = dir(fullfile(ekfPath, 'ekf*.data'));
            if ~isempty(ekfFilefindingDir)
                jointDataPath = fullfile(ekfPath, ekfFilefindingDir.name);
            else
                jointDataPath = [];
            end
        end
        
        % now for the end eff
        ekfFilefindingDir = dir(fullfile(ekfPath, 'EndEffector*.data'));
        if ~isempty(ekfFilefindingDir)
            endEffDataPath = fullfile(ekfPath, ekfFilefindingDir.name);
        else
            endEffDataPath = [];
        end
        
        % copy 
        if length(currExer) > 12
            currExerType = currExer(1:12);
        else
            currExerType = currExer(1:8);
        end

        cropPath = fullfile(currExerPath, 'Segmentation_cropping');
        cropHeaderPath = fullfile(cropPath, 'SegmentData_Cropping_Manual.header');
        cropDataPath = fullfile(cropPath, 'SegmentData_Cropping_Manual.data');
        
        manSegPath = fullfile(currExerPath, 'Segmentation_manual');
        if exist(manSegPath, 'file')
            manSegHeaderPath = fullfile(manSegPath, 'SegmentData_Manual_Manual.header');
            manSegDataPath = fullfile(manSegPath, 'SegmentData_Manual_Manual.data');
            
%             subjectInfo = [];
%             continue
        else
            % find other manual segments, perhaps they are renamed
            % differently
            manSegExtraDir = dir(fullfile(currExerPath, 'Segmentation_manual*'));
            
            if ~isempty(manSegExtraDir)
                manSegPath = fullfile(currExerPath, manSegExtraDir.name);
                manSegHeaderPath = fullfile(manSegPath, 'SegmentData_Manual_Manual.header');
                manSegDataPath = fullfile(manSegPath, 'SegmentData_Manual_Manual.data');
            else
                manSegHeaderPath = [];
                manSegDataPath = [];
                
%                 subjectInfo = [];
%                 continue
            end
        end
        
        % load device information
        [shimmerInfo] = shimmerLoadHeader(shimmerPath, shimmerDevicesFiles);
        
        segmentParam.date = datestr(now);
        segmentParam.segName = currExer;
        deviceInfo.deviceId{1} = shimmerInfo{1}.btId;
        deviceInfo.deviceId{2} = shimmerInfo{2}.btId;
        deviceInfo.deviceId{3} = shimmerInfo{3}.btId;
        deviceInfo.deviceLocation = shimmerDevicesFiles;
        
        % load shimmer data
        
        
        % load cortex data if exist
        if exist(cortexPath, 'file')
            cortexData = parseCSV_cortex(cortexPath);
        else
            cortexData = [];
        end
       
        % load joint angle data if exist
        if ~isempty(jointDataPath)
            % load joint angles data
            jointAngleData = parseCSV_jointAngle(jointDataPath);
            
            dt = mean(diff(jointAngleData.SystemSTimeStamp));
            
            % filter the data
            QTemp = jointAngleData.Q;            
            jointAngleData.Q = filter_dualpassBW(QTemp);
            dQTemp = diff(jointAngleData.Q)/dt;
            jointAngleData.dQ = dQTemp([1 1:size(dQTemp, 1)], :);
            ddQTemp = diff(jointAngleData.dQ)/dt;
            jointAngleData.ddQ = ddQTemp([1 1:size(ddQTemp, 1)], :);
        else
            jointAngleData.Q = [];
            jointAngleData.dQ = [];
            jointAngleData.ddQ = [];
        end
        
        if ~isempty(endEffDataPath)
            % load end effector data
            endEffectorData = parseCSV_endEffector(endEffDataPath);
            
            endEffectorData.EF1 = filter_dualpassBW(endEffectorData.EF1);
            endEffectorData.EF2 = filter_dualpassBW(endEffectorData.EF2);
        else
            endEffectorData.EF1 = [];
            endEffectorData.EF2 = [];
        end
       
        % load existing segmentation data
        if exist(cropDataPath, 'file')
            segmentData_crop = parseCSV_segment(cropDataPath);
        else
            segmentData_crop = [];
        end
        
        if exist(manSegDataPath, 'file')
            segmentData_manSeg = parseCSV_segment(manSegDataPath);
        else
            segmentData_manSeg = [];
        end
        
        % check the integrity of the data
        
        %if there is any nan or if it is too short, skip the data
        if sum(sum(isnan(jointAngleData.Q))) > 0 || size(jointAngleData.Q, 1) < 30 
            continue
        end
        
        subjectInfoInd = subjectInfoInd + 1;
        
        subjectInfo{subjectInfoInd}.subject = currSubjNum;
        subjectInfo{subjectInfoInd}.session = currSessNum;
        subjectInfo{subjectInfoInd}.exerciseName = currExer;
        subjectInfo{subjectInfoInd}.exerciseType = currExerType;
        subjectInfo{subjectInfoInd}.filePath = currExerPath;
        subjectInfo{subjectInfoInd}.time = jointAngleData.SystemSTimeStamp;
        subjectInfo{subjectInfoInd}.Q = jointAngleData.Q;
        subjectInfo{subjectInfoInd}.dQ = jointAngleData.dQ;
        subjectInfo{subjectInfoInd}.ddQ = jointAngleData.ddQ;
        subjectInfo{subjectInfoInd}.EF1 = endEffectorData.EF1;
        subjectInfo{subjectInfoInd}.EF2 = endEffectorData.EF2;
        subjectInfo{subjectInfoInd}.cortex = cortexData;
        subjectInfo{subjectInfoInd}.segmentData_crop = segmentData_crop;
        subjectInfo{subjectInfoInd}.segmentData_manSeg = segmentData_manSeg;
    end
end