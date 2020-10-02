function [obj, T, position] = calcOptInit(obj, pOptParam)
    % use FMINCON to minimize distance between observed marker position and
    % the forward kinematic position of the marker. 
    % x_opt = [tran_x, trans_y, trans_z, rot_x, rot_y, rot_z, jointparams]

    % parse the input struct and setup the optimization and initial est
    [optParam, x0_array] = parseParamAndInit(obj, pOptParam);
    
    % perform a fit with all the markers available
    optParam.markerWeight = ones(size(obj.model.sensors));
    [x_opt1, fval1] = runOpt(obj, x0_array, optParam);
    [minVal, minInd] = min(fval1);
 
    x_opt = x_opt1{minInd};
    [J, output1] = placement_func(obj, x_opt, optParam);
%     J_metric = max(cell2mat(output1(:, 2)));
    
    if 0
%         [J, output] = placement_func(obj, x_opt, markerWeight);
        [obj, T, position] = updateModel(obj, x_opt, optParam);
        plot_modelMarkers(obj);
        
        [obj, T, position] = updateModel(obj, x0_start, optParam);
        plot_modelMarkers(obj);
    end
    
    [obj, T, position] = updateModel(obj, x_opt, optParam);
end

function [optParam, x0_array] = parseParamAndInit(obj, pOptParam)
    p = inputParser;
    p.KeepUnmatched = true;

    addOptional(p, 'optTransXYZ', zeros(3, 1)); 
    addOptional(p, 'optRotXYZ', zeros(3, 1)); 
    addOptional(p, 'optJoints', zeros(length(obj.model.joints), 1)); 
    addOptional(p, 'frameInd', 1);
    addOptional(p, 'markersIndTransXYZ', []); 
    addOptional(p, 'markersIndRotXYZ', []); 

    parse(p, pOptParam);
    optParam = p.Results;
    
    % generate initial estimate for translational XYZ for the world-to-start frame
    if isempty(optParam.markersIndTransXYZ)
        initTransXYZ = [0 0 0]';
    else
        markerCum = [];
        for i = 1:length(optParam.markersIndTransXYZ)
            marker = obj.measurement_indArray{optParam.markersIndTransXYZ(i)}(optParam.frameInd, 1:3);
            markerCum = [markerCum; marker];
        end
        initTransXYZ = mean(markerCum, 1);
    end
    
    % generate initial estimate for rotation RPY for the world-to-start frame
    if isempty(optParam.markersIndRotXYZ)
        initRotXYZ = [0 0 0]';
    else
        % initialization for angle
        marker1 = obj.measurement_indArray{optParam.markersIndRotXYZ(1)}(optParam.frameInd, 1:2);
        marker2 = obj.measurement_indArray{optParam.markersIndRotXYZ(2)}(optParam.frameInd, 1:2);
        markers = [marker1; marker2];
        midPt = mean(markers);
        initAng = -atan2(midPt(2), midPt(1));
        
        initRotXYZ = [0 0 initAng]';
    end
    
    initJoints = zeros(length(obj.model.joints), 1);
    
    % now construct an init estimate for fmincon
    x0 = [];
    x0_char = {};
    for i = 1:length(optParam.optTransXYZ)
        if optParam.optTransXYZ(i) == 1
            x0 = [x0 initTransXYZ(i)];
            x0_char = [x0_char ['optTrans' num2str(i)]];
        end
    end
    
    for i = 1:length(optParam.optRotXYZ)
        if optParam.optRotXYZ(i) == 1
            x0 = [x0 initRotXYZ(i)];
            x0_char = [x0_char ['optRotat' num2str(i)]];
        end
    end
   
    for i = 1:length(optParam.optJoints)
        if optParam.optJoints(i) == 1
            x0 = [x0 initJoints(i)];
            x0_char = [x0_char ['optJoint' num2str(i)]];
        end
    end
    
    x0_array{1} = x0;
    
    % generate an array of rotArray if we're rotating the base yaw
    if optParam.optRotXYZ(3) == 1
        rotArray = pi/4:pi/4:2*pi;
        indToChange = find(strcmpi(x0_char, 'optRotat3'));
        x0_orig = x0;
        
        for i = 1:length(rotArray)
            x0 = x0_orig;
            x0(indToChange) = x0(indToChange) + rotArray(i);
            x0_array{end+1} = x0;
        end
    end
    
    optParam.x0_char = x0_char;
end

function [x_opt, fval] = runOpt(obj, x0_array, param)
    Mesoptions = optimset('Algorithm','interior-point','FunValCheck','on',...
        'MaxFunEvals',1e8,'MaxIter',100, 'Display','none');
  
    cost_func_handle = @(x_opt) placement_func(obj, x_opt, param);
    
    for i = 1:length(x0_array)
        x0 = x0_array{i};
        [x_opt{i}, fval(i)] = fmincon(cost_func_handle, x0, [],[],[],[],[],[],[],Mesoptions);
    end
end

function [J, output] = placement_func(obj, x0, optParam)
    % cost function
    [obj, T, position] = updateModel(obj, x0, optParam);
    error = zeros(length(obj.model.sensors), 1);
    
    for i = 1:length(obj.model.sensors)
        markerToUse = i;
        modelPos = obj.model.sensors(markerToUse).transform(1:3, 4);
        obsPos = obj.measurement_indArray{markerToUse}(1, 1:3)';
        baseError = norm(modelPos - obsPos, 2)*1000; % report the error in [cm] to inflate the sq fct
        error(i) = optParam.markerWeight(i)*baseError^2;
        
        output{i, 1} = obj.model.sensors(markerToUse).name;
        output{i, 2} = error(i);
    end

    J = sum(error);
end

function [obj, T, position] = updateModel(obj, x0, optParam)
    T = eye(4);
    position = zeros(length(obj.model.joints), 1);
    
    for i = 1:length(optParam.x0_char)
        operationString = optParam.x0_char{i}(1:8);
        operationCount = str2num(optParam.x0_char{i}(9:end));
        
        switch operationString
            case 'optTrans'
                T(operationCount, 4) = T(operationCount, 4) + x0(i);
                
            case 'optRotat'
                switch operationCount
                    case 1
                        T(1:3, 1:3) = rotx(x0(i))*T(1:3, 1:3);
                    case 2
                        T(1:3, 1:3) = roty(x0(i))*T(1:3, 1:3);
                    case 3
                        T(1:3, 1:3) = rotz(x0(i))*T(1:3, 1:3);
                end
                
            case 'optJoint'
                position(operationCount) = position(operationCount) + x0(i);
        end
    end
    
    obj.model.transforms(1).t = T;
    obj.model.position = wrapToPi(position);
    obj.model.forwardPosition();
end