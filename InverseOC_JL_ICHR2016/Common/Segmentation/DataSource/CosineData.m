classdef CosineData < AData
    % Cosine data generation
    
    % testing code
%     blah = CosineData(0:pi/32:8*pi, 1, 0);
%     blah.plot
%     blah.partition(5)
    
    properties(SetAccess=protected)
        dataNoiseless = [];
        data = []; % full array
        label = [];
        settings = [];
        
        trainingData = []; % subset, training purposes
        trainingLabel = [];
        
        testingData = []; % subset, testing purposes
        testingLabel = [];
        
        label_notSegment = 0;
        label_segment = 1;
    end
    
    properties(SetAccess=public)
        % public classes
        trainingDataDR = [];
        testingDataDR = [];
    end
    
    methods        
        function obj = CosineData(dataSettings,t,amp,phase)
            fprintf('Loading cosine dataset \n');
            
            obj.settings = dataSettings;
            obj.label_notSegment = obj.settings.label_notSegment;
            obj.label_segment = obj.settings.label_segment;
            
            obj.init(t,amp,phase);
        end
        
        function init(obj,t,amp,phase)
            % Initialization function 
            % Input:
            % t - time array
            % amplitude - 
            % phase - 
            
            % generate the cosine function
            obj.dataNoiseless = amp*cos(t - phase)';
            
            % add Gaussian noise
            noiseSNR = 25; % 1 is a lot of nose. 100 is no noise
            obj.data = awgn(obj.dataNoiseless, noiseSNR, 'measured');
            
            % generate labeling data based on zero crossings
            [crossingStruct, ~] = zvcCheck(t, obj.data');
            
            % expand the zero crossing windows so that it's '3' away from
            % each crossing point
            segmentPt = segmentPointWindowExpand(crossingStruct{1}.Index, obj.settings.segmentPointWindow);
            
            obj.label = ones(size(obj.data, 1), 1) * obj.label_notSegment;
            obj.label(segmentPt) = obj.label_segment;
        end
        
        function plot(obj)
            figure
            plot(obj.data)
            hold on
            crossingInd = find(obj.label == obj.label_segment);
            scatter(crossingInd, zeros(size(crossingInd)), 'rx');
        end
        
        function val = arraySize(obj)
            val = size(obj.data);
        end
        
        function partition(obj)
            % partition populates the training and testing arrays based on
            % window length
            
            % pull out all the 'segment' and 'not-segment' entries
            [trainingDataSegment, trainingLabelSegment] = pullDataForTrainingCosine(obj, obj.label_segment, obj.settings.dimStrack);
            [trainingDataNotSegment, trainingLabelNotSegment] = pullDataForTrainingCosine(obj, obj.label_notSegment, obj.settings.dimStrack);
            
            obj.trainingData = [trainingDataSegment; trainingDataNotSegment];
            obj.trainingLabel = [trainingLabelSegment; trainingLabelNotSegment];    
            
            % set the obs
            [testingDataAll, testingLabelAll] = pullDataForTestingCosine(obj, obj.settings.dimStrack);
            obj.testingData = testingDataAll;
            obj.testingLabel = testingLabelAll;
        end
    end
end