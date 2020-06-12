classdef BoostAggregator < AAggregator
% QDA Classifier class
    properties(SetAccess=protected)
       iterCount = 3;
       
       classifierEnsemble = [];
       classifierAlpha = [];
    end
    
    methods        
        function obj = BoostAggregator(classifier)
            obj.base_classifier = classifier;
        end
        
        function init(obj,iterCount)
            % Initialization function for Boosting algorithm
            obj.iterCount = iterCount;
        end
        
        function error = train(obj,input,output)
            % train the classifier
            classifier_array = cell(1, obj.iterCount);
            classifierLabel = cell(1, obj.iterCount);
            alpha = zeros(1, obj.iterCount);
            
            % initialize the weights. equal weights initially
            D = ones(length(output), 1)/length(output);
            
            indToUse = 1:length(output);
            
            iterD = D;
            iterInput = input;
            iterOutput = output; % these will be resampled per iteration
            
            Dtemp = [];
            
            for ind_iter = 1:obj.iterCount
                % train and apply classifier
                classifier_array{ind_iter} = obj.base_classifier.copy();   %Copy base
                classifier_array{ind_iter}.train(iterInput, iterOutput); %Train new
                classifierLabel{ind_iter} = classifier_array{ind_iter}.classify(iterInput); %Classify
                
                score = accuracyAssess(classifierLabel{ind_iter}, iterOutput);
                 
                % calculate classifier score, err
                err = sum(iterD .* score.incorrect_array);
                
                % calculate classifier weight, alpha
                alpha(ind_iter) = 0.5*log((1-err)/(err));
                
                if err == 0
                    % no error. current classifier is perfect
                    alpha(ind_iter) = 1;
                    break;
                end
                
                % rebalance training data weights, iterD
%                 iterD = iterD .* exp(-alpha(ind_iter) .* output .* classifierLabel{ind_iter});
                iterD = iterD .* exp(-alpha(ind_iter) .* score.incorrect_array);
                
                % update the D matrix with the results from iterD
                D(indToUse) = iterD;
                D = D ./ sum(D);
                
                % resample training data according to the weights
                indToUse = resampleBasedOnWeights2(D);
%                 indToUse = obj.resampleBasedOnWeights(D);
                
                if isempty(indToUse)
                    % some failsafe
                    indToUse = 1:length(output);
                end
                
                iterD = D(indToUse);
                iterInput = input(indToUse, :);
                iterOutput = output(indToUse, :); % these will be resampled per iteration
                
                Dtemp = [Dtemp; D'];
            end
            
            obj.classifierEnsemble = classifier_array;
            obj.classifierAlpha = alpha;
            
            error = [];
        end
        
        function out = classify(obj,input)
            H = zeros(size(input, 1), size(obj.classifierEnsemble, 2));
            for ind_iter = 1:size(obj.classifierEnsemble, 2)
                classifier = obj.classifierEnsemble{ind_iter};
                alpha = obj.classifierAlpha(ind_iter);
                
                if ~isempty(classifier)
                    H(:, ind_iter) = alpha * classifier.classify(input);
                end
            end
            
            out = abs(sign(sum(H, 2))); % not sure why it's coming out neg
        end
    end
    
    methods(Static,Hidden)
        function newDInd = resampleBasedOnWeights(D)
            % not sure how to do this. if it has a heavier weight,
            % include it more often, as par Seiffert2008

            % note that the returned newDInd is a function of the input D vector,
            % which may not be the full length of the original training data

            % pick 1/(n+1+wp) amounts of new points to keep, based on
            % the weights. change the wp as needed
            weightPenalty = 1;

            DUnique = unique(D); %
            DInd = [];

            % figure out how many entries are in each 'stage'
            DCount = zeros(size(DUnique));
            for ind_Dweights = 1:length(DUnique)
                DCount(ind_Dweights) = length(find(D == DUnique(ind_Dweights)));
            end

            % take the entry with the highest weights, and calculate resampling
            % ratios from that one
            countHeaviest = DCount(end);
            DUse = zeros(size(DUnique));
            for ind_Dweights = 1:length(DUnique)
                denom = length(DUnique)-ind_Dweights+1+weightPenalty;
                toUse = floor(countHeaviest/denom);

                if toUse > DCount(ind_Dweights)
                    % in case there is not enough entries in the lower weights
                    toUse = DCount(ind_Dweights);
                end

                DUse(ind_Dweights) = toUse;

                newD = find(D == DUnique(ind_Dweights));
                resampleInd = randperm(DUse(ind_Dweights)); % randomly pick some indices
                newDToUse = newD(resampleInd(1:DUse(ind_Dweights)));

        %         for ind_rep = 1:ind_Dweights
                    DInd = [DInd; newDToUse];
        %         end
            end

            newDInd = sort(DInd);
        end
        
    end
end