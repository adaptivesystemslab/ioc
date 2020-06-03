function basePath = dataBasePath(projectName)

    dataRoot = fullfile('D:', 'MyDocument', 'MotionData');

        switch projectName
            case 'Healthy1'
                projectPath = 'Lowerbody_healthy1_2011-11';

            case 'TRO'
                projectPath = '';

            case 'TRI'
                projectPath = 'Lowerbody_TRI1_2012-10';
        end

    basePath = fullfile(dataRoot, projectPath);

