function settings = pollDataSettings(name)
    % expects a sequence of strings, separated by underscores. 
    % param 1: dim stack - the width of window used for each timestep
    % param 2: exp - the number of p0 points to convert to p1 points around
    % the manual segment points
    % param 3: the feature to use. expects a string of 0s and 1s to
    % determine what is activated and what is not is determined in
    % dataGeneral (search 'FEATURES TO USE'). 
    
    % example: 15_10_110 will be n_stack = 15, n_exp = 10, with q and dq
    % enabled
    
    % dimStack - exp - q, dq, ddq, ef, def 15_10_110 15_10_000001
    
    parsing = regexp(name, '_', 'split');
    
    settings.dimStrack = str2num(parsing{1});
    settings.segmentPointWindow = str2num(parsing{2}); % how much non-segment points to convert to segment point
    
    switch parsing{3}
        case '0'
            % normal. load all data available
            settings.mansegLoad = 'loadAll';
            
%         otherwise
%              % semi-supervised. load the first segment, then try to learn 
%             settings.mansegLoad = ['loadFirst' parsing{3}];
    end
    
    settings.dataNormalization = parsing{4}; % zero the joint angles, normalize joint angles, abs the joint angles, normalize the joint velo
    
    jointVectorString = parsing{5};
    jointVectorString2 = regexp(jointVectorString, '-', 'split');
    
    if length(jointVectorString2) == 1
        % legacy mode
        jointVector = zeros(length(jointVectorString), 1);
        for ind = 1:length(jointVectorString)
            jointVector(ind) = str2num(jointVectorString(ind));
        end
    else
        for ind = 1:length(jointVectorString2)
            if ~isempty(jointVectorString2{ind})
                vectorInd = str2num(jointVectorString2{ind});
                jointVector(vectorInd) = 1;
            end
        end
    end
    
    settings.jointVector = jointVector;
    % convert it to string in case we're not on the legacy system
    jointVectorString = [];
    for ind = 1:length(jointVector)
        if jointVector(ind) == 1
            jointVectorString = [jointVectorString '1'];
        else
            jointVectorString = [jointVectorString '0'];
        end
    end
    
    settings.settingName = [parsing{1} '_' parsing{2} '_' parsing{3} '_' parsing{4} '_' jointVectorString];
    
     settings.segmentPointWindowVeloThreshold = 2.0;
end