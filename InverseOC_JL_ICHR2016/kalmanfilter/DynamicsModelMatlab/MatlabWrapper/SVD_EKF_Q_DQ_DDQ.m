classdef SVD_EKF_Q_DQ_DDQ < EKF_Q_DQ_DDQ
    %This is EKF that uses a Kinematic Model with sensors and can have
    %constraints on state and measurement using SVD
    
    properties     
       upper_state_bound = [];
       lower_state_bound = [];
       upper_mes_bound = [];
       lower_mes_bound = [];
    end
    
    methods
        function obj = SVD_EKF_Q_DQ_DDQ(model)
           obj = obj@EKF_Q_DQ_DDQ(model);
           
           %Init Bounds
           obj.upper_state_bound = ones(obj.sizeX,1).*inf;
           obj.lower_state_bound = ones(obj.sizeX,1).*-inf;
           
           obj.upper_mes_bound = ones(obj.sizeZ,1).*inf;
           obj.lower_mes_bound = ones(obj.sizeZ,1).*-inf;
        end
    end
    
    methods(Access=protected)
       
        function K = makeK(obj)
            %Call regular make K
            K = obj.makeK@EKF_Q_DQ_DDQ();
            K_new = K;
            %This is regular state update which we want bounded
            x = obj.x_predict+K*obj.dz;
            
            %---------Now apply bounds on the state using SVD-------------%
            [U,S,V] = svd(K);
            
            const_upper = x > obj.upper_state_bound;
            %This is q dqformulation of EKF, if q constraint is not
            %sattisfied impose constraint on dq to be 0
            dof = numel(obj.model_handle.joints);
            constraints = obj.upper_state_bound;
            I = find(const_upper(1:dof));
            if ~isempty(I)
                const_upper(I+dof) = x(I+dof) > 0;
                constraints(I+dof) = 1e-5;
            end
            
            %Project state prediction vector component by component onto the principle comps
            X_upper = U'*diag(obj.x_predict);
            %Project bounds component by component onto the principle comps
            B_mat_upper = U'*diag(constraints);            
            %Now go through each projection and see if the bounds are sattisfied
            y = V'*obj.dz;
            alphas = zeros(size(S,1),numel(x));
            
            for i=1:size(X_upper,2)
                if(const_upper(i))
                    %The bound is not sattisfied we need to adjust eigen values 
                    S_new = S;
                    
                   for j=1:size(S,1)
                      alpha = (B_mat_upper(j,i)-X_upper(j,i))./((S(j,:)*y));
                      if isnan(alpha)
                         1; 
                      end
                      if isinf(alpha) || isnan(alpha) 
                         alpha = 1; 
                      end
                      if alpha > 1
                         alpha = 1; 
                      end
                      if alpha < -1
                         alpha = -1; 
                      end

                      S_new(j,:) = S(j,:)*alpha;
                      alphas(j,i) = alpha;
                   end   
                   
                   K_new(i,:) = U(i,:)*S_new*V';
                   x=obj.x_predict+K_new*obj.dz;            %Update State Estimate
                   const_upper = x > constraints; %Update Constraints
                end
            end
            
            %Same idea for lower state bound
            const_lower = x < obj.lower_state_bound;
            
            %This is q dqformulation of EKF, if q constraint is not
            %sattisfied impose constraint on dq to be 0
            constraints = obj.lower_state_bound;
            I = find(const_upper(1:dof));
            if ~isempty(I)
                const_lower(I+dof) = x(I+dof) > 0;
                constraints(I+dof) = -1e-5;
            end
            
            
            X_lower = U'*diag(obj.x_predict);
            B_mat_lower = U'*diag(constraints);
            for i=1:size(X_lower,2)
                if(const_lower(i))
                    %The bound is not sattisfied we need to adjust eigen values 
                   S_new = S;
                    
                   for j=1:size(S,1)

                      alpha = (B_mat_lower(j,i)-X_lower(j,i))./((S(j,:)*y));
                      if isnan(alpha)
                         1; 
                      end
                      if isinf(alpha) || isnan(alpha) 
                         alpha = 1; 
                      end
                      if alpha > 1
                         alpha = 1; 
                      end
                      if alpha < -1
                         alpha = -1; 
                      end

                      S_new(j,:) = S(j,:)*alpha;
                      alphas(j,i) = alpha;
                   end    
                   
                   K_new(i,:) = U(i,:)*S_new*V';
                   x=obj.x_predict+K_new*obj.dz;            %Update State Estimate
                   const_lower = x < constraints; %Update Constraints
                end
            end
            K = K_new;
            
            %-------------------------STATE BOUNDS FINISHED--------------%
%             %-------------------------APPLY SIMILAR IDEA FOR PREDICTED
            %MEASUREMENTS------------------------------------------------%
            H_p = pinv(obj.H);
            const_upper = (obj.z_predict.getMesArray() - obj.upper_mes_bound) > 0;
            Y_upper = U'*H_p*diag((obj.upper_mes_bound-obj.z_predict));
            for i=1:size(Y_upper,2)
                if(const_upper(i))
                    %The bound is not sattisfied we need to adjust eigen values 
                   for j=1:size(S,1)

                      alpha = Y_upper(j,i)./((S(j,:)*y));
                      if isinf(alpha) || isnan(alpha) 
                         alpha = 0; 
                      end
                      if alpha > 1
                         alpha = 1; 
                      end
                      if alpha < -1
                         alpha = -1; 
                      end

                      S(j,:) = S(j,:)*alpha;
                      alphas(j,i) = alpha;
                   end       
                                       %Update 
                    K = U*S*V';
                    x=obj.x_predict+K*obj.dz;                           %Update State Estimate

                    obj.model_handle.position = x(1:dof);
                    obj.model_handle.velocity = x(dof+1:end);
                    obj.model_handle.forwardKinematics();

                     %H = obj.makeH(x);
                     z_new = obj.makeMeasure(x);

                     const_upper = (z_new - obj.upper_mes_bound) > 0;
                     Y_upper = U'*H_p*diag((obj.upper_mes_bound-obj.z_predict));
                     
                    %The bound is not sattisfied we need to adjust eigen values 
%                    for j=1:size(S,1)
%                       q = (Y_upper(j,i)./y(j))-S(j,j);
%                       S(j,j) = S(j,j) + q;
%                    end       
                end
            end
            
            const_lower = (z_new - obj.lower_mes_bound) < 0;
            Y_lower = U'*H_p*diag((obj.lower_mes_bound-obj.z_predict));
            for i=1:size(Y_lower,2)
                if(const_lower(i))
                    %The bound is not sattisfied we need to adjust eigen values 
                   for j=1:size(S,1)

                      alpha = Y_lower(j,i)./((S(j,:)*y));
                      if isinf(alpha) || isnan(alpha) 
                         alpha = 1; 
                      end
                      if alpha > 1
                         alpha = 1; 
                      end
                      if alpha < -1
                         alpha = -1; 
                      end

                      S(j,:) = S(j,:)*alpha;
                      alphas(j,i) = alpha;
                   end   
                   
                   %Update 
                    K = U*S*V';
                    x=obj.x_predict+K*obj.dz;                           %Update State Estimate

                    obj.model_handle.position = x(1:dof);
                    obj.model_handle.velocity = x(dof+1:end);
                    obj.model_handle.forwardKinematics();

                     %H = obj.makeH(x);
                     z_new = obj.makeMeasure(x);

                     const_lower = (z_new - obj.lower_mes_bound) < 0;
                     Y_lower = U'*H_p*diag((obj.lower_mes_bound-obj.z_predict));
                end
            end 
        end  
    end
end