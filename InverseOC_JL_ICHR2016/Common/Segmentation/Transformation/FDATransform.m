classdef FDATransform < ATransform
% Pricipal Component Analysis transform class

    properties (SetAccess=protected)
        %The different components
        coeffs = [];
        %How many components are used (default 1)
        degree = 1;
        elbow = [];
    end

    methods    
        function obj = FDATransform(degree)
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
        % Train FDA Transform

            %Get unique classes
            classes = unique(output);

            means = zeros(size(input,2),numel(classes));
            Sw = zeros(size(input,2),size(input,2));
            for i=1:numel(classes)
               means(:,i) = mean(input(output==classes(i),:)); 
               Sw = Sw + cov(input(output==classes(i),:)); 
            end

            %Find the total covariance
            St = cov(input);
            Sb = St-Sw;

            %Now find the eigenvectors
            [obj.coeffs, S, ~] = svd(Sw\Sb);
            
            %Calculate Elbow if we chose to
            if obj.degree == -1
                eig_vals = diag(S);
                eig_vals = eig_vals/sum(eig_vals);
                eig_diff = diff(eig_vals);
                sum_eig = 0;
                index = 0;
                
                %We maximum use 99% of the data or untill difference in
                %eigenvalues stops changing
                while sum_eig<0.99
                    sum_eig = sum_eig + eig_vals(index+1);
                    index = index +1;
                    
                    if(index > 4)
                       if(sum(eig_diff(index-4:index))<0.001)
                           obj.elbow = index-4;
                           break;
                       end
                    end
                    
                end
                obj.elbow = index;
            end
            
            error = 0;
        end
        
        function out = apply(obj,input)
            if(obj.degree == -1)
                out = (input*obj.coeffs(:,1:obj.elbow));
            else
                out = (input*obj.coeffs(:,1:obj.degree));
            end
        end
    end
end