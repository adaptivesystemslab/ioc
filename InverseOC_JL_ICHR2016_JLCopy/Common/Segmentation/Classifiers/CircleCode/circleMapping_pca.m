% a script package that loads data, generates PCA, generates a warping
% matrix for each subject/primitive pair, then calculates the inter-matrix
% warping distance
clearvars

loadExistingData = 1;

dataSelect.p0_guassianDownsampling = 1;
dataSelect.halfSegments = 1;
dataSelect.allowThresholdMultiplier = 1; % multiplier allowed on unbalanced datasets

dataSelect.outputBase = 'Results_Seg_SVMInc'; % root folder for output

classifierPackage.dataSelectList = {'15_10_1_1101_110'}; % dimStack - exp - q, dq, ddq, ef, def 15_10_110 15_10_000001
classifierPackage.dimReductionList = {'PCA_0'};
classifierPackage.aggregatorList = {'None'};
classifierPackage.classifierList = {'SVM_radial_0'};

batchSettings.batchInstancePath = 'None';
batchSettings.exportPathSuffix = 'None';

% load the training data
if ~loadExistingData
    % template training data
%     pxSet = 1;
%     count = 0;
%     sharcnetCommonSetupScript_script8_3;
%     trainingPackage = {};
%     testingPackage = {};
%     trainingPackage{1} = trainingPackageTemp; % we're not actually doing any training
%     
%     dataSelect.name = 'General';
%     trainingStruct = batchClassifier_2(batchSettings, dataSelect, trainingPackage, testingPackage, classifierPackage);
%     save trainingStruct
    
    pxSet = 1;
    count = 0;
    sharcnetCommonSetupScript_script8_1;
    testingPackage = {};
    trainingPackageLa{1} = trainingPackage(1); % we're not actually doing any training
    trainingPackageLa{2} = trainingPackage(5); % we're not actually doing any training
    trainingPackageLa{3} = trainingPackage(9); % we're not actually doing any training
    trainingPackageLa{4} = trainingPackage(13); % we're not actually doing any training
    trainingPackageLa{5} = trainingPackage(17); % we're not actually doing any training
    
    for i = 1:5
        dataSelect.name = 'General';
        trainingStruct_pca{i} = batchClassifier_2(batchSettings, dataSelect, trainingPackageLa{i}, testingPackage, classifierPackage);
    end
    save trainingStruct_pca
    
%     % testing data
%     pxSet = 1;
%     count = 0;
%     sharcnetCommonSetupScript_script8_3;
%     trainingPackage = {};
%     testingPackage = {};
%     trainingPackage{1} = testingPackageTemp; % we're not actually doing any training
%     
%     dataSelect.name = 'General';
%     testingStruct = batchClassifier_2(batchSettings, dataSelect, trainingPackage, testingPackage, classifierPackage);
%     save testingStruct
else
    load trainingStruct
%     load testingStruct
end

% abstract out the important variables
subjData = trainingStruct.trainingData.subjectData;
dimReduct = trainingStruct.dimReduct;

% parse through the training data and create a mapping matrix for each 
existingSubjInfo.subjectNumber = 0;
existingSubjInfo.sessionNumber = 0;
existingSubjInfo.exerciseName = '';
existingSubjInfo.exerciseType = '';

% this is the data for the templates
counter = 1;
for ind_subjData = 1:length(subjData)
    % is this the same subject as the last time?
    currSubj = subjData{ind_subjData};
    
    currSubjInfo.subjectNumber = currSubj.subjectNumber;
    currSubjInfo.sessionNumber = currSubj.sessionNumber;
    currSubjInfo.exerciseName = currSubj.exerciseName;
    currSubjInfo.exerciseType = currSubj.exerciseType;
    
    if strcmpi(currSubjInfo.exerciseName(end), '1')
        % do nothing
    else
        continue
    end
    
    counter = counter + 1;
    
    % extract each instance and apply stacking
    currData = [currSubj.q currSubj.dq];
    currLabel = currSubj.label;
    
    permissibleIndStruct = []; 
    rawTime = 1:length(currLabel);
    rawLabel = currLabel;
    rawSegTime = [];
    
    objSim.time = 1:length(currLabel);
    objSim.settings = pollDataSettings(classifierPackage.dataSelectList{1});
    objSim.settings.mode = 'Testing'; % need to create the proper I/Os
    objSim.label = currLabel;
    objSim.subjectData = subjData(ind_subjData);
    rawSegTestingInclude = currSubj.segTestingInclude;
    
    [testingData, testingLabel, testingInd, testingName, testingTime, testingSegTime, testingIndMapping] = ...
        samplingForTestingData(objSim, permissibleIndStruct, rawSegTime, [], rawSegTestingInclude);
    
    % then apply the proper PCA
    
    for ind_pca = 1:length(trainingStruct_pca)
        % find the proper pca
        if strcmpi(currSubj.exerciseName(1:8), trainingStruct_pca{i}.trainingData.motionSet{1}(1:8))
            dimReduct_pca = trainingStruct_pca{i}.dimReduct;
        end
    end
    pcaData = dimReduct_pca.apply(testingData); 
    
    % apply zero mean
%     pcaData = pcaData - repmat(mean(pcaData), size(pcaData, 1), 1); % remove the mean
    
    % fit a ellipse to the pca data to find the proper 0,0
%     ellipsePartStruct = fit_ellipse(pcaData(:, 1), pcaData(:, 2));
%     pcaData(:, 1) = pcaData(:, 1) - ellipsePartStruct.X0;
%     pcaData(:, 2) = pcaData(:, 2) - ellipsePartStruct.Y0;
%     
    % then convert to matrix
    [xp{counter}, P{counter}, h] = mapToCircle(pcaData, testingLabel, 1); % xp = Ppca' * (x-mu) * Pc;
    
    header{counter} = [currSubj.datasetName '~' ...
        'Subj' num2str(currSubj.subjectNumber) '~' ...
        'Sess' num2str(currSubj.sessionNumber) '~' ...
        currSubj.exerciseName];
    sortHeader{counter} = [currSubj.exerciseName '~' ...
        'Subj' num2str(currSubj.subjectNumber)];
    
    if ~isempty(h)
        title(header{counter});
        saveas(h, ['C:\Documents\MATLABResults\imgDump\circlemapping\nominal_' sortHeader{counter} '.jpg']);
        close(h);
    end
end

[sortMag, sortInd] = sort(sortHeader);
header = header(sortInd);
P = P(sortInd);
table = cell(size(header, 2)+1, size(header, 2)+1);
table(2:end, 1) = header';
table(1, 2:end) = header;

% calculate distance metric 
for ind_i = 1:size(P, 2)
    for ind_j = 1:size(P, 2)
        % lets calc us some norms
        
        A = P{ind_i}(:, 1:2);
        B = P{ind_j}(:, 1:2);
        
        d_1 = matrixNorm(A, B, 1);
        d_2 = matrixNorm(A, B, 2);
        d_inf = matrixNorm(A, B, Inf);
        d_cos = matrixNorm(A, B, 'cos');
        d_fro = matrixNorm(A, B, 'fro'); % identical to matlab's norm(A, 'fro')
        
        d_matlab_1 = norm(A - B, 1);
        d_matlab_2 = norm(A - B, 2);
        d_matlab_inf = norm(A - B, Inf);
        d_matlab_cos = 0;
        d_matlab_fro = norm(A - B, 'fro');
        
        comparison = [d_1 d_2 d_inf d_cos d_fro; 
            d_matlab_1 d_matlab_2 d_matlab_inf d_matlab_cos d_matlab_fro];
        
        mtxDist_1norm(ind_i, ind_j) = d_1;
        mtxDist_2norm(ind_i, ind_j) = d_2;
        mtxDist_inf(ind_i, ind_j) = d_inf;
        mtxDist_cos(ind_i, ind_j) = d_cos;
        mtxDist_fro(ind_i, ind_j) = d_fro;
    end
end

mat2cellVal = size(mtxDist_1norm);
table_1norm = table; table_1norm(2:end, 2:end) = mat2cell(mtxDist_1norm, ones(mat2cellVal(1), 1), ones(1, mat2cellVal(2)));
table_2norm = table; table_2norm(2:end, 2:end) = mat2cell(mtxDist_2norm, ones(mat2cellVal(1), 1), ones(1, mat2cellVal(2)));
table_inf = table;   table_inf(2:end, 2:end) = mat2cell(mtxDist_inf, ones(mat2cellVal(1), 1), ones(1, mat2cellVal(2)));
table_cos = table;   table_cos(2:end, 2:end) = mat2cell(mtxDist_cos, ones(mat2cellVal(1), 1), ones(1, mat2cellVal(2)));
table_fro = table;   table_fro(2:end, 2:end) = mat2cell(mtxDist_fro, ones(mat2cellVal(1), 1), ones(1, mat2cellVal(2)));