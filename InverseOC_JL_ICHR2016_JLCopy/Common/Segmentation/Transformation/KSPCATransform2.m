classdef KSPCATransform2 < ATransform
    % ksPCA, reparsed from Ghodsi
    
%  Usage: k = kernel(ker,u,v,p1,p2)
%
%  Parameters: ker - kernel type
%              u,v - kernel arguments
%
%  Values for ker: 'linear'     - 
%                  'delta'      -  
%                  'poly'       - p1 is degree of polynomial
%                  'rbf'        - p1 is width of rbfs (sigma)
%                  'sigmoid'    - p1 is scale, p2 is offset
%                  'spline'     -
%                  'bspline'    - p1 is degree of bspline
%                  'fourier'    - p1 is degree
%                  'erfb'       - p1 is width of rbfs (sigma)
%                  'anova'      - p1 is max order of terms
    
    properties(SetAccess=private)
        
       %kernel functions
       y_kernal = 'none';
       x_kernal = 'none';
       
       y_kernal_param1 = [];
       y_kernal_param2 = [];
       x_kernal_param1 = [];
       x_kernal_param2 = [];
       
       trainingData = [];
       
       sigma = [];
        
       mean_x = [];
       mean_y = [];
       
       U = [];
       D = [];
       alignRot = eye(3);
       wmean = [0 0 0]';
       
       degree = 2;
       LMatrixMultiplier = 1;
       
       downsampleUpperBound = 2500; % so n p1 points, and n p0 points
    end
        
    methods
        function obj = KSPCATransform2()
           obj.y_kernal = 'delta';
           obj.x_kernal = 'none'; 
        end
        
        function init(obj,varargin)
           %Initialization input:
           %    Kernel function y: 'linear', 'rbf', 'delta'
           %    Kernel function x: 'linear', 'rbf', 'none'
           %    sigma
           %    degree
           
           obj.y_kernal = varargin{1};
           obj.x_kernal = varargin{2};
           obj.sigma = varargin{3};
           obj.degree = varargin{4};    
           obj.LMatrixMultiplier = varargin{5};
        end
        
        function error = train(obj,X,Y)
        
            % apply downsample and save the training data
            if size(X, 1)/2 > obj.downsampleUpperBound
                % randomly select some points for downsampling
                p1points = find(Y == 1); % don't want to assume the array is sorted
                p0points = find(Y == 0);
                
                selectInd = randperm(length(p1points), obj.downsampleUpperBound); % first half is p1 points
                selectInd2 = randperm(length(p0points), obj.downsampleUpperBound);
                
                X_tr = X([p1points(selectInd) p0points(selectInd2)], :)';
                Y_tr = Y([p1points(selectInd) p0points(selectInd2)], :)';
            else
                X_tr = X';
                Y_tr = Y';
            end
            
            % the paper requires data to be normalized
%             for ind_length = 1:(size(X_tr, 1))
%                 holding = X_tr(ind_length, :) - min(X_tr(ind_length, :));
%                 X_tr(ind_length, :) = holding / max(holding);
%             end
                
%             [Y_tr, Y_tr_sortInd] = sort(Y_tr, 'descend'); % sort the array so can use heuristics to generate kernel matrices
%             X_tr = X_tr(:, Y_tr_sortInd);

            obj.trainingData = X_tr; % save the array for future usage
            
            param.LMatrixMultiplier = obj.LMatrixMultiplier;
            param.ktype_y = obj.y_kernal;
            param.ktype_x = obj.x_kernal;
            
            if obj.degree == -1
                % run PCA to get the cutoff
                testPCA = PCATransform(-1);
                testPCA.train(X_tr',Y_tr');
                obj.degree = testPCA.elbow;
            end
             
            param.kparam_y = 1; % 1 in the example code
            d = obj.degree;
                    
            % set up parameters
            switch obj.x_kernal
                case 'none'
                   
                case 'rbf'
                    param.kparam_x = obj.sigma;

                otherwise
                    param.kparam_x = obj.degree;
            end
            
            % call sPCA code
            switch obj.x_kernal
                case 'none'
%                     Y_tr(Y_tr == 1) = 2;
%                     Y_tr(Y_tr == 0) = 1;
                    [Ztr_SPCA U D] = SPCA(X_tr,Y_tr,d,param);
                    
                otherwise
                    [Ztr_KSPCA U D] = KSPCA(X_tr,Y_tr,d,param);
            end
            
            obj.U = U;
            obj.D = D;
            
            % 4 mvt
            
            % 5 mvt
            
                
%             % actually, use this other instead 
% healthy1 5 mvt (set1, KHEF in other direction)

% healthy1 5mvt (set2, KHEF in the same direction)
% scaled version with eye*1 of set 2
     
% healthy1 kefo only
% 
%             obj.U = Uprime;
            
            
            
%             obj.alignRot = eye(size(U, 2));
%             obj.wmean = zeros(size(U, 2), 1);
            
%             cylindicalAlignment(obj, obj.apply(X), Y);
        
        end
        
        function cylindicalAlignment(obj, tData, tLabel)

         
         % calculate regression onto a flat plane
         tData_x = tData(:, 1);
         tData_y = tData(:, 2);
         tData_z = tData(:, 3);
         
         pointLimit = 5000;
         points_p1 = findpPoints(tLabel, 1, pointLimit); % all the points that are p1, then randomly selected for plotting
         points_p0 = findpPoints(tLabel, 0, pointLimit); % all the points that are p0, then randomly selected for plotting

         if 0
          h = figure;             hold on
tData_p1 = tData(points_p1, :);
            tData_p0 = tData(points_p0, :);
            
            scatter3(tData_p1(:, 1), tData_p1(:, 2), tData_p1(:, 3), ...
                'CData', [0 0 1])

            % plot p0 points
            scatter3(tData_p0(:, 1), tData_p0(:, 2), tData_p0(:, 3), ...
                'CData', [1 0 0])

        grid on; view([-20 30]); 
            xlim([-15 15]); ylim([-15 15]); zlim([-15 15])
            xlabel('x'); ylabel('y'); zlabel('z')

         end
         
         left(1, 1) = sum(tData_x .^ 2);
         left(1, 2) = sum(tData_x .* tData_y);
         left(1, 3) = sum(tData_x);
         left(2, 1) = left(1, 2);
         left(2, 2) = sum(tData_y .^ 2);
         left(2, 3) = sum(tData_y);
         left(3, 1) = left(1, 3);
         left(3, 2) = left(2, 3);
         left(3, 3) = length(tData_x);
         
         right(1) = sum(tData_x .* tData_z);
         right(2) = sum(tData_y .* tData_z);
         right(3) = sum(tData_z);
         
         soln = left \ right'; 
         
         % project all datapoints onto calculated plane
%          [X,Y] = meshgrid(-10:10:10, -10:10:10);
%          Z = soln(1)*X + soln(2)*Y + soln(3);
%          surf(X,Y,Z);
% xlabel('x'); ylabel('y');
         
%          normVector = [soln(1)/soln(3), soln(2)/soln(3), -1/soln(3)];
%          v = [tData_x tData_y tData_z - soln(3)];
         
        normVector = [soln(1), soln(2) -1];
        v = [tData_x tData_y tData_z];
      
         for ind = 1:length(v)
             vind = v(ind, :);
             w(ind, :) = vind - normVector*dot(vind, normVector); % w is the projected points
         end
         
         % shift centre of the w so that it is at [0, 0, 0]
         wmean = mean(w);
         soln = soln - wmean';
         w = w-repmat(wmean, size(w, 1), 1);
         

         % now rotate the plane so it's flat in cartesian, so it's easier
         % to work with
         
         % rotate in y
         X = normVector(1);
         Z = normVector(3);
         angToRot2 = atan(X/Z);
         Ry = roty(angToRot2);
         blarg = Ry'*normVector';
        
         X = blarg(1);
         Y = blarg(2);
         Z = blarg(3);
         angToRot2 = atan(Y/Z);
         Rx = rotx(-angToRot2);
         blarg2 = Rx'*blarg;

         
fulldata = ((Ry*Rx)'*(v-repmat(wmean, size(w, 1), 1))')';
         
%          Y = 0;
%          X = -10;
%          Z = soln(1)*X + soln(2)*Y + soln(3);
%          angToRot2 = atan(Z/X);
%          Ry = roty(angToRot2);
%          ack = (Ry'*w')';
%          soln = Ry'*soln;
         
           % rotate in x
%          X = 0;
%          Y = 10;
%          Z = soln(1)*X + soln(2)*Y + soln(3);
%          angToRot3 = atan(Z/Y);
%          Rx = rotx(angToRot3);
%          ack = (Rx*ack')';
%          soln = Rx*soln;
         
Rz = eye(3);

% % % % % % pull out all the p1 points 
% % % % % fulldata_p1points = fulldata(points_p1, :);
% % % % %     [cidx, ctrs] = kmeans(fulldata_p1points, 2, 'Distance','city');
% % % % % %                           figure; grid on; hold on
% % % % % %            plot3(fulldata_p1points(cidx==1,1),fulldata_p1points(cidx==1,2),fulldata_p1points(cidx==1,3),'r.', ...
% % % % % %              fulldata_p1points(cidx==2,1),fulldata_p1points(cidx==2,2),fulldata_p1points(cidx==2,3),'b.', ctrs(:,1),ctrs(:,2),ctrs(:,3),'kx');
% % % % %                           
% % % % %          % now rotate in z
% % % % %          X = ctrs(1);
% % % % %          Y = ctrs(2);
% % % % %          angToRot3 = atan(Y/X);
% % % % %          Rz = rotz(angToRot3);
% % % % % %          ack = (Rz*ack')';
% % % % %          soln = Rz*Ry*soln;
         
         obj.alignRot(1:3, 1:3) = Ry*Rx*Rz;
         obj.wmean(1:3) = wmean';
         
         
         polb = (obj.alignRot(1:3, 1:3))'*normVector';
         
         %          figure; hold on; grid on
%          xlabel('x'); ylabel('y'); zlabel('z');
%          plot3([0 normVector(1)], [0 normVector(2)], [0 normVector(3)], 'b')
         %          plot3([0 blarg(1)], [0 blarg(2)], [0 blarg(3)], 'r')
 %          plot3([0 blarg2(1)], [0 blarg2(2)], [0 blarg2(3)], 'g')
%          plot3([0 polb(1)], [0 polb(2)], [0 polb(3)], 'k')

         fulldata2 = ((obj.alignRot(1:3, 1:3))'*(v-repmat(wmean(1:3), size(w, 1), 1))')';
         
%          scatter3(v(:, 1), v(:, 2), v(:, 3))
%          
%          figure; hold on; grid on;
%          xlabel('x'); ylabel('y'); zlabel('z');
% %          scatter3(fulldata1(:, 1), fulldata1(:, 2), fulldata1(:, 3))
%          scatter3(fulldata2(:, 1), fulldata2(:, 2), fulldata2(:, 3))
% 
%             tData_p1 = tData(points_p1, :);
%             tData_p0 = tData(points_p0, :); hold on
%             scatter3(tData_p1(:, 1), tData_p1(:, 2), tData_p1(:, 3), 'b')
%             scatter3(tData_p0(:, 1), tData_p0(:, 2), tData_p0(:, 3), 'r'); view([20 50])
% 
%    tData_p1 = fulldata2(points_p1, :);
%             tData_p0 = fulldata2(points_p0, :); hold on
%             scatter3(tData_p1(:, 1), tData_p1(:, 2), tData_p1(:, 3), 'b')
%             scatter3(tData_p0(:, 1), tData_p0(:, 2), tData_p0(:, 3), 'r')


        end
        
        function out = apply(obj,input)
            
            X_tr = obj.trainingData;
            X_ts = input';
            
%             % normalizing obs data
%             for ind_length = 1:(size(X_ts, 1))
%                 holding = X_ts(ind_length, :) - min(X_ts(ind_length, :));
%                 X_ts(ind_length, :) = holding / max(holding);
%             end
            
            param.ktype_y = obj.y_kernal;
            param.ktype_x = obj.x_kernal;
             
            param.kparam_y = 1; % 1 in the example code
            d = obj.degree;
            
            switch obj.x_kernal
                case 'none'
                    out = obj.U'*X_ts;
                    
                otherwise
                    switch obj.x_kernal
                        case 'rbf'
                            param.kparam_x = obj.sigma;
                        otherwise
                            param.kparam_x = obj.degree;
                    end
                    
                    Ktest = zeros(size(X_tr,2),size(X_ts,2));
                    for i=1:size(X_tr,2)
                        for j=1:size(X_ts,2)
                            Ktest(i,j) = kernel(param.ktype_x,X_tr(:,i),X_ts(:,j),param.kparam_x,[]);
                        end
                    end
                    out = obj.U'*Ktest;
            end

            out = out';
            
%             out = ((obj.alignRot)'*(out-repmat(obj.wmean', size(out, 1), 1))')';
            
            
        end
    end
        
    
    
end