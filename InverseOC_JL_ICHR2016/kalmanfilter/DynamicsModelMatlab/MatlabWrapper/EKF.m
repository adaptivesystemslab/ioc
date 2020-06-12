classdef EKF < handle
    %This is EKF implemented with the very basic equations
    
    properties
        observation_noise = [];
        process_noise = [];
        covariance = [];
        state = [];
    end
    
    properties (SetAccess = protected)
        isWriting = false;
        %Size of state
        sizeX = 0;
        %size of Measurement
        sizeZ = 0;
        
        %Jacobians
        H = [];
        A = [];
        
        %Measurement Error Covariance
        inov_covariance = [];
        inov_covariance_inv = [];
        %Predicted State Covariance
        P_predict;
        
        %Predicted State
        x_predict;
        
        %Predicted Measurement
        z_predict;
        
        %Measurement Residual
        dz;
        
        %Optimal Kalman Gain
        K = [];
        
        % number of iterations performed
        iterationStep = 0;
    end
    
    methods
        
        function init(obj,file_path)
            %Initializes the file writer and creates the header file for the
            %EKF.
            %@param: path full path of the file if file already exists it will
            %be cleared
            if(exist([file_path '.header'],'file') || exist([file_path '.data'],'file'))
                %error('EKF: File already exists');
                delete([file_path '.header']);
                delete([file_path '.data']);
            end
            
            %Write XML header
            docNode = com.mathworks.xml.XMLUtils.createDocument('EKF');
            docRootNode = docNode.getDocumentElement;
            docRootNode.setAttribute('date',datestr(clock));
            
            %Process Noise
            el = docNode.createElement('ProcessNoise');
            str = sprintf('%f, ',diag(obj.process_noise));
            el.appendChild(docNode.createTextNode(str(1:end-2)));
            docRootNode.appendChild(el);
            
            %Observation Noise
            el = docNode.createElement('ObservationNoise');
            str = sprintf('%f, ',diag(obj.observation_noise));
            el.appendChild(docNode.createTextNode(str(1:end-2)));
            docRootNode.appendChild(el);
            
            xmlFileName = [file_path,'.header'];
            xmlwrite(xmlFileName,docNode);
            
            obj.data_file = fopen([file_path '.data'],'w');
            
            %Print the headers
            header_str = 'SystemMSTimeStamp,Dt,';
            for i=1:obj.sizeX
                header_str = strcat(header_str,['S' num2str(i) ',']);
            end
            fprintf(obj.data_file,'%s\r\n',header_str);
            obj.isWriting = true;
            
        end
        
        function run_iteration(obj,u,z,timestamp)
            %Run EKF Iteration
            
            obj.iterationStep = obj.iterationStep + 1;
            
            obj.x_predict = obj.stateUpdate(obj.state,u);              %State Update Equation
            obj.A = obj.makeA(obj.x_predict,u);                        %State Update Jacobian
            obj.H = obj.makeH(obj.x_predict);                          %Measurement Jacobian
            
            obj.P_predict = obj.A*obj.covariance*obj.A'+obj.process_noise;%partial update
            obj.inov_covariance = obj.H*obj.P_predict*obj.H' + obj.MakeObservationNoise();              %cross covariance
            
            obj.z_predict = obj.makeMeasure(obj.x_predict);        %Measurement Prediction
            obj.dz = obj.makeDZ(z,obj.z_predict);                       %Measurement Residual
            
            obj.K = obj.makeK();
            
            obj.state=obj.x_predict+obj.K*obj.dz;                           %Update State Estimate
            obj.covariance = obj.makeP();                           %Updated Covariance Estimate
            
            %obj.write(u,timestamp);
        end
                
        function z = makeMeasure(obj,x)
            %Function to generate Measurement from State should be
            %overridden by any implementation, default assumes state is
            %measured
            
            z = x;
            
        end
        
        function dz = makeDZ(obj,z_act,z_predict)
            %Measurement residual, default is simple difference
            dz = z_act-z_predict;
        end
        
        function H = makeH(obj,x)
            %Make the observation Jacobian, should be overridden in any EKF implementation
            %The default generates the Jacobian numerically
            
            %Get predicted measurement given current state
            z = obj.makeMeasure(x);
            z_numeric = false;
            if isa(z,'numeric')
                H = zeros(numel(z),numel(x));
                z_numeric = true;
            else
                H = zeros(numel(z.getMesArray),numel(x));
            end
            %We modify the state a bit at a time and get new measurement
            eps = 1e-8;
            for i=1:numel(x)
                %Change state a little bit
                x(i) = x(i) + eps;
                %Get new measurement
                z_new = obj.makeMeasure(x);
                if(z_numeric)
                   H(:,i) = (z_new-z)/eps;
                else
                   H(:,i) = (z_new.getMesArray-z.getMesArray)/eps; 
                end
                %for j=1:numel(z)
                %    H(j,i) = (z_new(j)-z(j))/eps;
                %end
                %Change state back
                x(i) = x(i)-eps;
            end
        end
        
        function x_new = stateUpdate(obj,x,u)
            x_new = x;
        end
       
        
        function A = makeA(obj,x,u)
            %Generate state update Jacobian, should be overridder in derived
            %classes, default is the identity matrix
            
            A = eye(numel(x));
            
        end
        
        function mes_cov = projectCov(obj)
            %Projected measurement covariance
            mes_cov = obj.H*obj.covariance*obj.H';
        end
        
    end
    
    methods(Access=protected)
        %Protected Methods of EKF
        
        function K = makeK(obj)
            %Make Optimal Kalman Gain
            K = obj.P_predict*obj.H'/obj.inov_covariance;            %Optimal Kalman Gain
        end
        
        function x = makeState(obj)
            x=obj.x_predict+obj.K*obj.dz;                       %Update State Estimate
        end
        
        function value = MakeObservationNoise(obj)
            value = obj.observation_noise;
        end
        
        function P = makeP(obj)
            %Update covariance using Kalman Gain, we use the Joseph Form to
            %maintain PSD
            
            %Regular Form
            %P = obj.P_predict-obj.K*obj.H*obj.P_predict;
            
            %Joseph Form
            try
            I_m_KH = eye(size(obj.K,1)) - obj.K*obj.H;
            P = I_m_KH*obj.P_predict*I_m_KH' + obj.K*obj.MakeObservationNoise()*obj.K';
            catch
               disp('asdf') 
            end
        end
        
        function write(obj,dt,timestamp)
            if(obj.isWriting == true)
                if nargin > 3
                    fprintf(obj.data_file,'%f,',timestamp);
                else
                    fprintf(obj.data_file,'%f,',java.lang.System.currentTimeMillis);
                end
                fprintf(obj.data_file,'%f,',dt);
                fprintf(obj.data_file,'%f,',obj.state);
                fprintf(obj.data_file,'%s\r\n','');
            end
        end
        
    end
    
end