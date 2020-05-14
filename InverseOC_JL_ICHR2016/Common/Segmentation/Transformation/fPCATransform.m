classdef fPCATransform < ATransform
% Pricipal Component Analysis transform class

    properties (SetAccess=protected)
        %The different components
        PCAModel = [];
        
        threshold = 0;
    end

    methods    
        function obj = fPCATransform(degree, threshold)
        % Constructor for PCA Transform
        % Input: degree
        %       Number of principal components data will be projected on
            if ~exist('threshold', 'var')
               obj.threshold = 0.8;
            else
                obj.threshold = threshold;
            end
            
%             obj.PCAModel = PCATransform(3, obj.threshold);
        end
        
        function init(obj,varargin)
        % Initialization of PCA Transform
        % Input: degree
        %       Number of principal components data will be projected on
            obj.degree = varargin(1);
        end
        
        function error = train(obj,dataStruct)
        % use the real PCA to determine the number of degrees to use
        tempPCA = PCATransform(-1, obj.threshold);   
        tempPCA.train(dataStruct.trainingData, dataStruct.trainingLabel);        
        obj.PCAModel = PCATransform(tempPCA.elbow, tempPCA.threshold);
        
        % convert the fake PCA data into something usable
        fakeDataFrame = dataStruct.dataFake;
        dt = mean(diff(dataStruct.time(1:10)));
        t = 0:dt:5-dt;
        fakeData = [];
        
        for i = 1:5
            % pull the max/min behaviour of the data observed
            maxDofVal = max(fakeDataFrame(:, i));
            minDofVal = min(fakeDataFrame(:, i));
            
            tx1 = Septic(minDofVal,0,0,0,maxDofVal,0,0,0,t);
            tx2 = Septic(maxDofVal,0,0,0,minDofVal,0,0,0,t);
            fakeData = [fakeData [tx1 tx2]'];
        end
        
        fakeData = [fakeData [0 0 0 0 0; diff(fakeData)/dt]];
        
%         for i = 6:10
%             % pull the max/min behaviour of the data observed
%             maxDofVal = max(fakeDataFrame(:, i))*2;
%             minDofVal = min(fakeDataFrame(:, i))*2;
%             
%             tx1 = Septic(minDofVal,0,0,0,maxDofVal,0,0,0,t(1:length(t)/2));
%             tx2 = Septic(maxDofVal,0,0,0,minDofVal,0,0,0,t(1:length(t)/2));
%             fakeData = [fakeData [tx1 tx2 -tx2 -tx1]'];
%         end
        
        genTime = [t t+t(end)+dt]';
        genInput = fakeData;
        genLabel = ones(size(genTime)); % we're training PCA so label shouldn't matter
        
        % perform the dim stacking        
        subjData{1}.subjectNumber = 0;
        subjData{1}.data = genInput;
        subjData{1}.label = genLabel;
        subjData{1}.name = [];
        
        objSim.time = 1:length(genLabel);
        objSim.settings = dataStruct.settings;
        objSim.settings.mode = 'Testing'; % need to create the proper I/Os
        objSim.label = genLabel;
        objSim.subjectData = subjData;
        
        [testingData, testingLabel, testingInd, testingName, testingTime, testingSegTime, testingIndMapping] = ...
            samplingForTestingData(objSim, [], [], [], genLabel);
        
        % call PCA
        
        
        obj.PCAModel.train(testingData, testingLabel);
        end
        
        function out = apply(obj,input)
%             input = bsxfun(@minus, input, obj.mapping.mean);
%             out = input * obj.mapping.M;
            out = obj.PCAModel.apply(input);
        end
    end
end