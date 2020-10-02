function ekf = setupEKF(model, dt, ekfType, tuningParams)

    switch ekfType
        case 'MC_EKF_Q_DQ_DDQ'
            ekf = eval([ekfType '()']);
            ekf.eta = tuningParams.processNoiseCoefficient;
            for i = 1:length(model)
                ekf.addModel(model(i));
            end
          
            modelJoints = [model(1).joints model(2).joints];
            modelSensors = [model(1).sensors model(2).sensors];
            
        otherwise
            ekf = eval([ekfType '(model)']);   
            modelJoints = model.joints;
            modelSensors = model.sensors;
    end
    
    dim = numel(modelJoints);
    
    allSensTypes = join({modelSensors.type}, '');
    allTypes = strsplit(allSensTypes{1}, ',');
    allTypes = allTypes(1:end-1);
    allTypesUnique = unique(allTypes);
    
    obs_noise = [];
    for i=1:numel(modelSensors)
        sens = modelSensors(i);
        types = strsplit(sens.type,',');
        types = types(1:end-1);
        
        for j=1:numel(types)
            switch types{j}
                case 'position'
                    if sum(ismember(allTypesUnique, 'accelerometer'))
                        obs_noise = [obs_noise; ones(3,1)*...
                            tuningParams.observationNoiseCoefficientVirtualPosition];
                    else
                        obs_noise = [obs_noise; ones(3,1)*...
                            tuningParams.observationNoiseCoefficientMarkerPosition];
                    end
                    
                case 'velocity'
                    if sum(ismember(allTypesUnique, 'accelerometer'))
                        obs_noise = [obs_noise; ones(3,1)*...
                            tuningParams.observationNoiseCoefficientVirtualVelocity];
                    else
                        obs_noise = [obs_noise; ones(3,1)*...
                            tuningParams.observationNoiseCoefficientMarkerVelocity];
                    end
                    
                case 'accelerometer'
                    obs_noise = [obs_noise; ones(3,1)*...
                        tuningParams.observationNoiseCoefficientAccelerometer];
                    
                case 'gyroscope'
                    obs_noise = [obs_noise; ones(3,1)*...
                        tuningParams.observationNoiseCoefficientGyroscope];
            
                case 'yaw'
                    obs_noise = [obs_noise; ...
                        tuningParams.observationNoiseCoefficientYaw];
         
                otherwise
                    error(['Unknown Sensor Type: ' types{j}]);
            end
        end
    end
    
    obs_noise = diag(obs_noise);
    ekf.observation_noise = obs_noise;
    
    if isfield(tuningParams, 'processNoiseCoefficientPrismatic')
        etaPrismatic = tuningParams.processNoiseCoefficientPrismatic;
        etaRevolute = tuningParams.processNoiseCoefficientRevolute;
    else
        etaPrismatic = tuningParams.processNoiseCoefficient;
        etaRevolute = tuningParams.processNoiseCoefficient;
    end
    
    G = [];
    for i=ekf.sizeX/dim:-1:1
        G = [G; ones(dim,1)*(dt)^i / factorial(i)];
    end
    
    %apply eta sqrt
    stateLen = ekf.sizeX/dim;
    for i = 1:dim
        inds = [];
        for j = 1:stateLen
            inds = [inds i+(j-1)*dim]; 
        end
        
        switch model.joints(i).type
            case 'prismatic'
                G(inds) = G(inds) * sqrt(etaPrismatic);
                
            case 'revolute'
                G(inds) = G(inds) * sqrt(etaRevolute);
        end
    end
    
   %Process noise of EKF
    P_tmp = G*G';
    P = zeros(size(P_tmp));
    for i=1:ekf.sizeX/dim
        for j=i:ekf.sizeX/dim
            P(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim) = ...
                diag(diag(P_tmp(i*dim-dim+1:i*dim,j*dim-dim+1:j*dim)));
        end
    end

    P = P+P' - diag(diag(P));
    %P(1:dim,1:dim) = P(1:dim,1:dim) + eye(dim)*tuningParams.processNoiseCoefficientPosition;
    %P(dim+1:2*dim,dim+1:2*dim) =  P(dim+1:2*dim,dim+1:2*dim) + eye(dim)*tuningParams.processNoiseCoefficientVelocity;
    switch ekfType
        case 'MC_EKF_Q_DQ_DDQ'
            %Noise is set up for each model when the model is added
        otherwise
            ekf.process_noise = P;
    end
    
    ekf.covariance = tuningParams.covarianceInitialization*eye(ekf.sizeX);
end