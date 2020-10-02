classdef MC_EKF_Resize < EKF_Resize
 
    
    properties
        %List of model handles, we can add models to the EKF which will
        %update the size of state, cov, and mes, and ect
        %model_handle  = rlCModel.empty();
        
        %Keeps track of where the state starts for each model
        model_state_indices = 1;
        
        %Constraints on the models
        constraints = EKF_Constraint.empty();
    end
    
    methods
        
        
    end
    
    
end