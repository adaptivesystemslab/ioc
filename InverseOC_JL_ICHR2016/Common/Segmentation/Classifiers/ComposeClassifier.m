classdef ComposeClassifier < AClassifier
   %This is the interface function to the compose algorithm
   
    properties(Access=private)
       composeModel = [];
       modelParam = [];
       classifierName = [];
       
       trainingData = [];
       trainingLabel = [];
       
       ignoreVal = -Inf;
    end
    
    methods
        function obj = ComposeClassifier(varargin)
            switch str2num(varargin{1}.classifierName)
                case 1
                    obj.classifierName = 's3vm';
                case 2
                    obj.classifierName = 'label_prop';
                case 3
                    obj.classifierName = 'label_spread';
                case 4
                    obj.classifierName = 'cluster_n_label';
            end
        end
        
        function init(obj)
            
        end
        
        function error = train(obj,input, output)
           
            obj.trainingData = input;
            obj.trainingLabel = output;
        end
        
        function out = classify(obj,input,labels,incMode,extra)
            
            ind = 1;
            
            dataSet1 = zeroScaleSignal(obj.trainingData, 1);
%             dataSet1 = obj.trainingData;
            
            dataset.data{ind} = dataSet1(:, 1:3);
            dataset.labels{ind} = obj.remapLabels(obj.trainingLabel);
            dataset.use{ind} = ones(size(obj.trainingLabel)) * 1;
            dataset.labels{ind}(end) = 0;
            dataset.use{ind}(end) = 0;
            
            ind = 2;
            
            dataSet2 = zeroScaleSignal(input, 1);
%             dataSet2 = input;
            
            dataset.data{ind} = dataSet2(:, 1:3);
            dataset.labels{ind} = zeros(size(labels));
            dataset.use{ind} = ones(size(labels)) * 0;
            
            dataset.data = dataset.data';
            
            verbose = 2;
            
            obj.composeModel = compose(dataset, verbose);
            obj.composeModel.set_cores(4);
            obj.composeModel.set_classifier(obj.classifierName); %  {'s3vm'; 'label_prop'; 'label_spread'; 'cluster_n_label'};
            obj.composeModel.set_cse('a_shape');
            
            a = obj.composeModel.run;
            
            out = obj.revertLabels(obj.composeModel.hypothesis{2});
            
            acc = sum(out == labels)/length(labels);
            
        end
        
        function trainFsm(obj, trainingData)
            
        end
        
        function initializeAdaptiveElements(obj)
            
        end

        function obj_copy = copy(obj)
           obj_copy = SVMClassifier();
           obj_copy.init; % linear would = 0 if the +1 didn't exist
        end
    end    
    
    methods(Static)
        function labels = remapLabels( labels)
            labels(labels == 0) = 2;
        end
        
        function labels = revertLabels(labels)
            labels(labels == 2) = 0;
        end
    end
end