function [model, q, q_name, t] = loadPRMDATFiles(currSubjSDIMSFileName, currSubjAmnFileName)
    % given the path to the SDIMS file and the AMN folderpath, load the
    % joint angles
    
    % load the corresponding joint angles
    model = loadSDIMS(currSubjSDIMSFileName);
    jointAngles = [];
    jointAngleNames = {};
    
    for j = 3:length(model.Node)
        currNodeName = model.Node{j}.name;
        currPath = fullfile(currSubjAmnFileName, [currNodeName '.dat']); 
        currNodeQR = dlmread(currPath);
        t = currNodeQR(:, 1);
        d = currNodeQR(:, 2:end);
        
        switch size(d, 2)
            case 1
                % 1 DOF joint angle
                jointAngles = [jointAngles d];
                jointAngleNames{end+1} = [currNodeName '_q1'];
                
            case 9
                % 3 DOF rotational matrix
                jointAngles = [jointAngles reverseRotationalMatrix(d)];
                jointAngleNames{end+1} = [currNodeName '_q1'];
                jointAngleNames{end+1} = [currNodeName '_q2'];
                jointAngleNames{end+1} = [currNodeName '_q3'];
        end
    end
    
        % load the corresponding joint angles     
    model = loadSDIMS(currSubjSDIMSFileName);
    jointAngles = [];
    jointAngleNames = {};
    
    for j = 3:length(model.Node)
        currNodeName = model.Node{j}.name;
        currPath = fullfile(currSubjAmnFileName, [currNodeName '.dat']);
        currNodeQR = dlmread(currPath);
        t = currNodeQR(:, 1);
        d = currNodeQR(:, 2:end);
        
        switch size(d, 2)
            case 1
                % 1 DOF joint angle
                jointAngles = [jointAngles d];
                jointAngleNames{end+1} = [currNodeName '_q1'];
                
            case 9
                % 3 DOF rotational matrix
                jointAngles = [jointAngles reverseRotationalMatrix(d)];
                jointAngleNames{end+1} = [currNodeName '_q1'];
                jointAngleNames{end+1} = [currNodeName '_q2'];
                jointAngleNames{end+1} = [currNodeName '_q3'];
        end
    end
    
    q_name = jointAngleNames;
    q = jointAngles;
    
%     q_name = jointAngleNames(21:end);
%     q = jointAngles(:, 21:end);
    
%     q_filt = filter_dualpassBW(outAngles);