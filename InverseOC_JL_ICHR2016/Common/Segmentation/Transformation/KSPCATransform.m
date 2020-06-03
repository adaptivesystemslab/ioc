classdef KSPCATransform < ATransform
    
    properties(SetAccess=private)
        
       %kernel functions
       y_kernel = [];
       x_kernel = [];
       sigma = 0.5;
       
       downsampleUpperBound = 1000;
        
       mean_x = [];
       mean_y = [];
       
       U = [];
       training_points;
       degree = 2;
    end
        
    methods
        function obj = KSPCATransform()
           obj.y_kernel = 'delta';
           obj.x_kernel = 'none'; 
        end
        
        function init(obj,varargin)
           %Initialization input:
           %    Kernel function y: 'linear', 'rbf', 'delta'
           %    Kernel function x: 'linear', 'rbf', 'none'
           %    sigma
           %    degree
           
           obj.y_kernel = varargin{1};
           obj.x_kernel = varargin{2};
           obj.sigma = varargin{3};
           obj.degree = varargin{4};
           
        end
        
        function error = train(obj,input,output)
        
            if size(output, 1) > obj.downsampleUpperBound
                % randomly select some points for downsampling
                halfXLength = floor(size(input, 1)/2);
                selectInd = randperm(halfXLength, obj.downsampleUpperBound); % first half is p1 points
                selectInd2 = randperm(halfXLength, obj.downsampleUpperBound) + halfXLength;
                input = input([selectInd selectInd2], :);
                output = output([selectInd selectInd2], :);
            else
                input = input;
                output = output;
            end
            
            N = size(input,1);
            obj.mean_x = mean(input,1);
            obj.mean_y = mean(output,1);
            input = input-repmat(obj.mean_x,N,1);
            output = output-repmat(obj.mean_y,N,1);
            obj.training_points = input;
            %Now that we centered out matrix we want 
            %to find U which maximizes Tr(U'XBX'U)
            %Since we are kernalizing we actually assume U=X*Bt
            %So max Tr(Bt'X'XBX'XBt) notice that X'X = K is da kernel stuff
            
            B = obj.kernelmatrix(obj.y_kernel,output',[],obj.sigma);
            
            if(strcmp(obj.x_kernel,'none'))
                [u s v] = svd(input'*B*input,0);
            else
                K = obj.kernelmatrix(obj.x_kernel,input',[],obj.sigma);
                [u s v] = svd(K*B*K,0);
            end
            
            obj.U = u;
            
            %Calculate Elbow if we chose to
            if obj.degree == -1
                eig_vals = diag(s);
                eig_vals = eig_vals/sum(eig_vals);
                
                sum_eig = 0;
                index = 0;
                while sum_eig<0.80
                    sum_eig = sum_eig + eig_vals(index+1);
                    index = index +1;
                end
                obj.degree = index;
            end
        end
        
        function out = apply(obj,input)
            N = size(input,1);
            input = input-repmat(obj.mean_x,N,1);
            
            if(strcmp(obj.x_kernel,'none'))
                out = (obj.U(:,1:obj.degree)'*input')';
            else
                K = obj.kernelmatrix(obj.x_kernel,obj.training_points',input',obj.sigma);
                out = (obj.U(:,1:obj.degree)'*K)';
            end
        end
    end
        
    
    methods(Static)
        
        function K = kernelmatrix(ker,X,X2,sigma)
            switch ker
                case 'lin'
                    if exist('X2','var') && ~isempty(X2)
                        K = X' * X2;
                    else
                        K = X' * X;
                    end

                case 'poly'
                    if exist('X2','var') && ~isempty(X2)
                        K = (X' * X2 + b).^d;
                    else
                        K = (X' * X + b).^d;
                    end

                case 'rbf'
                    n1sq = sum(X.^2,1);
                    n1 = size(X,2);

                    if isempty(X2);
                        D = (ones(n1,1)*n1sq)' + ones(n1,1)*n1sq -2*X'*X;
                    else
                        n2sq = sum(X2.^2,1);
                        n2 = size(X2,2);
                        D = (ones(n2,1)*n1sq)' + ones(n1,1)*n2sq -2*X'*X2;
                    end;
                    K = exp(-D/(2*sigma^2));
                    
                case 'rbf_svkernel'
                    % X is either training data (for training data) or 
                    % X is training data + X2 is obs data
                    n1sq = sum(X.^2,1);
                    n1 = size(X,2);

                    if isempty(X2);
                        D = (ones(n1,1)*n1sq)' + ones(n1,1)*n1sq -2*X'*X;
                    else
                        n2sq = sum(X2.^2,1);
                        n2 = size(X2,2);
                        D = (ones(n2,1)*n1sq)' + ones(n1,1)*n2sq -2*X'*X2;
                    end;
                    K = exp(-D/(2*sigma^2));
                    
for i = 1:n
    for j = 1:n
        L(i,j) = L(i,j) + kernel(param.ktype_y,Y(:,i),Y(:,j),param.kparam_y,[]);
    end
end

                case 'sam'
                    if exist('X2','var') && ~isempty(X2)
                        D = X'*X2;
                    else
                        D = X'*X;
                    end
                    K = exp(-acos(D).^2/(2*sigma^2));
                    
                case 'delta'
                    N = size(X,2);
                    K = eye(N,N);

                    classes = unique(X)';

                    for i=1:size(classes,1)
                       K(X == classes(i,:),X == classes(i,:)) = 1; 
                    end
                otherwise
                    error(['Unsupported kernel ' ker])
            end
        end
       %Here are the kernel functions 
       function out = linear(input)
          out = input*input'; 
       end
       
       function out = radial(input)
           sigma = 0.1;
           input = input';
           n1sq = sum(input.^2,1);
           n1 = size(input,2);
           D = (ones(n1,1)*n1sq)' + ones(n1,1)*n1sq -2*input'*input;
           out = exp(-D/(2*sigma^2));
       end
       
       function out = delta(input)
          %Delta kernel for labels as in ghodsi's paper
          N = size(input,1);
          out = eye(N,N);
          
          classes = unique(input);
          
          for i=1:size(classes,1)
             out(input == classes(i,:),input == classes(i,:)) = 1; 
          end
       end
    end
    
end