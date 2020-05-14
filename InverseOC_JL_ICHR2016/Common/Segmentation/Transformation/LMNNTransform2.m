classdef LMNNTransform2 < ATransform
% LMNN using Weinberger's implementation
% http://www.cs.cornell.edu/~kilian/code/lmnn/lmnn.html

    properties (SetAccess=protected)
        %The different components
        U = [];
        m = [];
        %How many components are used (-1 means calculate the elbow of eigenvalues)
        degree = -1;
        elbow = [];
        
        L = [];
        Details = [];
        
%         dtwMtx = [];
%         dtwDeg = [];

        pcaModel = [];
        
        downsampleUpperBound = 50000; % so n p1 points, and n p0 points
    end

    methods    
        function obj = LMNNTransform2(degree)
        % Constructor for PCA Transform
        % Input: degree
        %       Number of principal components data will be projected on
            obj.degree = degree;
        end
        
        function init(obj,varargin)
        % Initialization of PCA Transform
        % Input: degree
        %       Number of principal components data will be projected on
            obj.degree = varargin(1);
        end
        
        function error = train(obj,input,output)

            % apply downsample and save the training data
            if size(input, 1)/2 > obj.downsampleUpperBound
                % randomly select some points for downsampling
                p1points = find(output == 1); % don't want to assume the array is sorted
                p0points = find(output == 0);
                
                selectInd = randperm(length(p1points), obj.downsampleUpperBound); % first half is p1 points
                selectInd2 = randperm(length(p0points), obj.downsampleUpperBound);
                
                X_tr = input([p1points(selectInd) p0points(selectInd2)], :);
                Y_tr = output([p1points(selectInd) p0points(selectInd2)], :)';
            else
                X_tr = input;
                Y_tr = output;
            end


        if obj.degree == 1
            obj.pcaModel = PCATransform(-1, 0.95);
            obj.pcaModel.train(X_tr,Y_tr);
            x = obj.pcaModel.apply(X_tr)';
        else
            x = X_tr';
        end
        
            y = Y_tr';
            
            numData = size(x, 2);
            numTraining = floor(numData/2);
            indTraining = sort(randperm(numData,numTraining));
            indValid = setxor(1:numData, indTraining);
            
            xTr=x(:,indTraining);
            xVa=x(:,indValid);
            yTr=y(:,indTraining);
            yVa=y(:,indValid);
            
        %% tune parameters
%         disp('Getting started');
        % These parameters were found with the following search command:
        [Klmnn,knn,outdim,maxiter]=findLMNNparams(xTr,yTr,xVa,yVa);
%         K=15;
%         outdim=68;
%         maxiter=120;
        
        %% train full muodel
%         fprintf('Training final model...\n');
        [obj.L, obj.Details] = lmnnCG([xTr xVa], [yTr yVa],Klmnn,'maxiter',maxiter,'outdim',outdim);
%         [testerr, testDetails] = knncl(obj.L,xTr,yTr,x,y,3,'train',0);
%         testerr=knncl(obj.L,xTr,yTr,xTe,yTe,3,'train',0);
%         fprintf('\n\nTesting error: %2.2f%%\n',100.*testerr);


if 0 
             h = figure;             hold on
          % plot p1 points
          x = x';
          tData_p1 = x(output == 1, :);
          tData_p0 = x(output == 0, :);

          title('PCA')
          scatter3(tData_p1(:, 1), tData_p1(:, 2), tData_p1(:, 3),  'CData', [0 0 1], 'Marker', '.')
          scatter3(tData_p0(:, 1), tData_p0(:, 2), tData_p0(:, 3), 'CData', [1 0 0], 'Marker', '.')

           x = obj.pcaModel.apply(input)';
            D=length(obj.L);
            x=obj.L*x(1:D,:);
            
            input2 = x';
            
            
                      h = figure;             hold on
          % plot p1 points
          title('LMNN')
          tData_p1 = input2(output == 1, :);
          tData_p0 = input2(output == 0, :);

          scatter3(tData_p1(:, 1), tData_p1(:, 2), tData_p1(:, 3),  'CData', [0 0 1], 'Marker', '.')
          scatter3(tData_p0(:, 1), tData_p0(:, 2), tData_p0(:, 3), 'CData', [1 0 0], 'Marker', '.') 
end
        end
        
        function out = apply(obj,input)
            if obj.degree == 1
                x = obj.pcaModel.apply(input)';
            else
                x = input';
            end
            
            D=length(obj.L);
            x=obj.L*x(1:D,:);
            
            out = x';
        end
        
        function out = counterApply(obj,input)
            x = input';
            
            D=length(obj.L);
            x=(obj.L'*x(1:D,:))';
            
            if obj.degree == 1
                x = obj.pcaModel.counterApply(x);
            else
                x = input;
            end
            
            out = x;
        end
    end
end