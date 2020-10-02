classdef Camera < matlab.mixin.Heterogeneous & handle

    % Protected C++ pointer
    properties (SetAccess = protected)
        % This is the actual integer pointer to C++ object in mem
        c_ptr           = 0;
        
        sensor;
        
        observation     = []; 
        projection      = [];
        window          = [];
        fullProjection  = [];
        
        features;
%         featureVector   = uint8(zeros([100, 64]));
        corners;
        count           = 0;
    end
    
    properties (Dependent = true)
        pose            = eye(4);
        calibration     = [671.7691,0,320;0,671.7691,180;0,0,1];
    end
    
    methods
        %Constructor  
        function obj = Camera(name, model, frame, preTransform)
            obj.c_ptr = DynamicsModelMatlab(11, 1, name);
            obj.sensor = SensorCore(name);

            if nargin < 4
                preTransform = eye(4);
            end
            
            model.addSensor(obj.sensor, frame, preTransform);
            model.forwardPosition;
            obj.updatePose;
        end
        
        function [features, corners] = getFeatures(obj)
            features = obj.features;
            corners = obj.corners;
        end
        
        function updatePose(obj)
            obj.pose = obj.sensor.transform;
        end
        
        function p = get.pose(obj)
            p = DynamicsModelMatlab(11, 4, obj.c_ptr);
            obj.pose = p;
        end
        
        function set.pose(obj, p)
            DynamicsModelMatlab(11, 3, obj.c_ptr, p);
        end
        
        function result = get.observation(obj)
            obj.observation = obj.getObservation(true);
            result = obj.observation;
        end
        
        function K = get.calibration(obj)
            K = DynamicsModelMatlab(11, 6, obj.c_ptr);
            obj.calibration = K;
        end
        
        function set.calibration(obj, K)
            d = size(K);
            if d(1) ~= 3 || d(2) ~= 3
                ME = MException('Camera:calibrationNot3x3', ...
                                'Camera Calibration Matrix needs to be 3x3.');
                throw(ME)
            end
            DynamicsModelMatlab(11, 5, obj.c_ptr, K);
        end
        
        function P = get.projection(obj)
            obj.projection = obj.getProjectionMatrix;
            P = obj.projection;
        end
        
        function W = get.window(obj)
            obj.window = obj.getWindowMatrix;
            W = obj.window;
        end
        
        function WP = get.fullProjection(obj)
            obj.fullProjection = obj.window * obj.projection;
            WP = double(obj.fullProjection);
        end
        
        function xw = inverseProjectPoint(obj, x_img, depth)
            T_0c = obj.pose;
            
            % Make Projection Matrix Invertible
            P = [obj.fullProjection; 0 0 -1];
            
            % Inverse Z Matrix.
            iZf = double(diag([-depth, -depth, -depth, 1]));
            
            % Invert the projection, then the homogenous coordinate
            % transformation, and finally go into world coordinates.
            x_pc = P \ [x_img; 1];           
            x_c = iZf * [x_pc; 1];
            xw = T_0c * x_c;
        end
        
        function xp = projectPoint(obj, x_world)
            T_0c = obj.pose;
            x_c = T_0c \ x_world;
            
            % (Big) Z matrix (represents the homogenous transform)
            Zf = [diag(repmat(-1/x_c(3), 1, 3)), [0;0;0;]];
            
            xp = obj.fullProjection * Zf * x_c;
        end
        
        function dxp = projectPointVelocity(obj, x_0, dq_dt)
            T_0c = obj.pose;
            T_c0 = [[T_0c(1:3, 1:3)', -T_0c(1:3, 1:3)'*T_0c(1:3,4)]; ...
                    [0, 0, 0, 1]];
            x_c = T_c0 * x_0;
            
            dT_c0 = obj.sensor.dTb_ee;
            dT_dt = zeros(4, 4);
            
            for i = 1:length(dq_dt)
                dT_dt = dT_dt + dq_dt(i) * dT_c0(:,:,i);
            end
            
            Zf = [diag(repmat(-1/x_c(3), 1, 3)), [0;0;0;]];
            dZf_dzc =[diag(repmat(1/((x_c(3))^2), 1, 3)), [0;0;0;]];
            K = obj.fullProjection;
            
            dT = dZf_dzc * ( dT_dt(3, :) * x_0 ) + Zf * dT_dt;
            dxp = K * dT * x_0;
        end
               
        function J = jacobian_q(obj, x_0)
            T_0c = obj.pose;
            T_c0 = [[T_0c(1:3, 1:3)', -T_0c(1:3, 1:3)'*T_0c(1:3,4)]; ...
                    [0, 0, 0, 1]];
            x_c = T_c0 * x_0;

            % Camera Projection Ma
            K = obj.fullProjection;
            
            % (Big) Z matrix (represents the homogenous transform)
            Zf = [diag(repmat(-1/x_c(3), 1, 3)), [0;0;0;]];
            dZf_dzc = [diag(repmat(1/((x_c(3))^2), 1, 3)), [0;0;0;]];
            
            dTc0_dq = obj.sensor.dTb_ee;
            
            % Joints in model is size of 3rd dim
            M = size(dTc0_dq,3);
            dXc = zeros(4, M);
            A = zeros(4, M);
            for i = 1:M
                dXc(:,i) = dTc0_dq(:,:,i) * x_0;
                A(:,i) = dXc(3,i) * x_c;
            end
            J = K * (Zf * dXc + dZf_dzc * A);
        end
        
        function J = jacobian_x0(obj, x_0)
            T_0c = obj.pose;
            T_c0 = [[T_0c(1:3, 1:3)', -T_0c(1:3, 1:3)'*T_0c(1:3,4)]; ...
                    [0, 0, 0, 1]];
            x_c = T_c0 * x_0;

            % Camera Projection Matrix
            K = obj.fullProjection;
            
            % (Big) Z matrix (represents the homogenous transform)
            Zf = [diag(repmat(-1/x_c(3), 1, 3)), [0;0;0;]];
            dZf_dzc = [diag([1/((x_c(3))^2),1/((x_c(3))^2),1/((x_c(3))^2)]),[0;0;0;]];
            
            J = [ K * (Zf * T_c0(:,1) + dZf_dzc * T_c0(3,1) * x_c), ...
                  K * (Zf * T_c0(:,2) + dZf_dzc * T_c0(3,2) * x_c), ...
                  K * (Zf * T_c0(:,3) + dZf_dzc * T_c0(3,3) * x_c) ];
              
%             % The Full Jacobian expanded out.
%             J = [ K * (Zf * T_c0 * [1;0;0;0] + dZf_dzc * ([0,0,1,0] * T_c0 * [1;0;0;0]) * T_c0 * x_0), ...
%                   K * (Zf * T_c0 * [0;1;0;0] + dZf_dzc * ([0,0,1,0] * T_c0 * [0;1;0;0]) * T_c0 * x_0), ...
%                   K * (Zf * T_c0 * [0;0;1;0] + dZf_dzc * ([0,0,1,0] * T_c0 * [0;0;1;0]) * T_c0 * x_0) ];
        end
        
        function result = getObservation(obj, waitForNewImage)
            % Retrieves the observed image. It will block until a fresh
            % image is available by default, so make sure to call vis.update 
            % before this if you do not pass an argument! 
            %
            % As a consequence of retrieving the image, the camera
            % will also cache the feature vector created by performing a
            % FAST analysis.

            if nargin < 2
                waitForNewImage = 1;
            end

            result = DynamicsModelMatlab(11, 7, obj.c_ptr, waitForNewImage);
            result = flip(result, 1);

            % Extract Corners from the Image and cache them.
            Ig = rgb2gray(result);
            specialPoints = detectMinEigenFeatures(Ig);
            specialPoints = specialPoints.selectStrongest(10);
            [obj.features, obj.corners] = extractFeatures(Ig, specialPoints);
        end

        function mat = getProjectionMatrix(obj)
            % This function returns the camera's projection matrix. The
            % Projection matrix transforms the frustum view cone of the camera
            % into a rectangular subspace. It is used to multiply a point in
            % camera space.
        
            % We need the calibration matrix to find out what the width and
            % height of the image is.
            K = obj.calibration;
            
            % Reduced Term Projection Matrix.
        	mat = diag([K(1,1)/K(1,3), -K(2,2)/K(2,3), -1]);

%             % The full OpenGL Matrix. We don't need it..and removed most
%             % components.
%             % The default near and far planes for the MRPT::CCamera class.
%             zn = 0.1;
%             zf = 10000;
%             
%             % Projection Matrix
%             mat = diag([K(1,1)/K(1,3), -K(2,2)/K(2,3), -(zf + zn)/(zf - zn), 0]);
%             mat(4,3) = -1;
%             mat(3,4) = -2*zf*zn/(zf-zn);
        end
        
        function mat = getWindowMatrix(obj)
            % The Window Matrix takes the normalized device coordinate point
            % and transforms it into pixel space. The normalized device coord
            % point is constructed by dividing the XYZ values by the 4th
            % component of the projection matrix output. That is:
            % 
            % Xc = obj.getProjectionMatrix * X
            % Xndc = [(Xc(1:2) / Xc(3); 1]
            % Xpixel = obj.getWindowMatrix * Xndc
            % 
            % More generally, the Window Matrix is a 4x4 Homogenous
            % Transform that simply scales the point and adds an offset.
        
            K = obj.calibration;
           
            mat = zeros([2, 3]);
            mat(1, 1) = K(1, 3);
            mat(2, 2) = -K(2, 3); % Flipped Y Axis
            mat(1:2,3) = K(1:2,3) + 1; % 1-Based Indexing for MATLAB
        end

        %Destructor
        function delete(obj)
            if obj.c_ptr ~= 0
                DynamicsModelMatlab(11, 2, obj.c_ptr);
            end
            obj.c_ptr = 0;
        end
   end
end