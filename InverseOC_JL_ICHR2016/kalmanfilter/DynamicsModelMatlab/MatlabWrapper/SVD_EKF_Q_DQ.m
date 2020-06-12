classdef SVD_EKF_Q_DQ < EKF_Q_DQ
    %This is EKF that uses a Kinematic Model with sensors and can have
    %constraints on state and measurement using SVD
    
    properties     
       upper_state_bound = [];
       lower_state_bound = [];
       upper_mes_bound = [];
       lower_mes_bound = [];
    end
    
    methods
        function obj = SVD_EKF_Q_DQ(model)
           obj = obj@EKF_Q_DQ(model);
           
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
            K = obj.makeK@EKF_Q_DQ();
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
                const_upper(I+dof) = 1;
                constraints(I+dof) = -1e-5;
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

                      
                      if isinf(alpha) || isnan(alpha) 
                         alpha = 0; 
                      end
%                       if alpha > 1
%                          alpha = 1; 
%                       end
%                       if alpha < -1
%                          alpha = -1; 
%                       end
                      S_new(j,:) = S(j,:)*alpha;
                      alphas(j,i) = alpha;
                   end   
                   
                   K_new(i,:) = U(i,:)*S_new*V';
                   %K_new = U*S_new*V';
                   x(i)=obj.x_predict(i)+U(i,:)*S_new*V'*obj.dz;            %Update State Estimate
                   const_upper = const_upper | (x > constraints); %Update Constraints
                end
            end
            
            %Update Partial 
            cols = const_upper == 0;
            H_tmp = obj.H(obj.selectedcols,cols);
            P_tmp = obj.P_predict(cols,cols);
            obs_noise = obj.MakeObservationNoise();
            inov_cov_tmp = H_tmp*P_tmp*H_tmp'+obs_noise(obj.selectedcols,obj.selectedcols);
            K_new(cols,:) = P_tmp*H_tmp'/inov_cov_tmp;
            
            %Same idea for lower state bound
            const_lower = x < obj.lower_state_bound;
            
            %This is q dqformulation of EKF, if q constraint is not
            %sattisfied impose constraint on dq to be 0
            constraints = obj.lower_state_bound;
            I = find(const_lower(1:dof));
            if ~isempty(I)
                const_lower(I+dof) = x(I+dof) > 0;
                constraints(I+dof) = 1e-5;
            end
            
            
            X_lower = U'*diag(obj.x_predict);
            B_mat_lower = U'*diag(constraints);
            for i=1:size(X_lower,2)
                if(const_lower(i))
                    %The bound is not sattisfied we need to adjust eigen values 
                    S_new = S;
                    
                   for j=1:size(S,1)
                      alpha = (B_mat_lower(j,i)-X_lower(j,i))./((S(j,:)*y));

                      if isinf(alpha) || isnan(alpha) 
                         alpha = 0; 
                      end
%                       if alpha > 100
%                          alpha = 100; 
%                       end
%                       if alpha < -100
%                          alpha = -100; 
%                       end
                      S_new(j,:) = S(j,:)*alpha;
                      alphas(j,i) = alpha;
                   end   
                   
                   K_new(i,:) = U(i,:)*S_new*V';
                   %K_new = U*S_new*V';
                   x=obj.x_predict+K_new*obj.dz;            %Update State Estimate
                   const_lower = const_lower | (x < constraints); %Update Constraints
                end
            end
            cols = const_lower == 0 & const_upper == 0;
            H_tmp = obj.H(obj.selectedcols,cols);
            P_tmp = obj.P_predict(cols,cols);
            obs_noise = obj.MakeObservationNoise();
            inov_cov_tmp = H_tmp*P_tmp*H_tmp'+obs_noise(obj.selectedcols,obj.selectedcols);
            K_new(cols,:) = P_tmp*H_tmp'/inov_cov_tmp;
            K = K_new;
            
            %-------------------------STATE BOUNDS FINISHED--------------%
            
%             %Run FK to get measurements 
%             dof = obj.model_handle.dof;
%             obj.model_handle.position = x(1:dof);
%             obj.model_handle.velocity = x(dof+1:end);
%             obj.model_handle.forwardKinematics();
%             
%              %H = obj.makeH(x);
%             z_new = obj.makeMeasure(x);

            %-------------------------APPLY SIMILAR IDEA FOR PREDICTED
            %MEASUREMENTS------------------------------------------------%
            
%             W = [1 1 1 0 0 0];
%             
%             %Pull out only columns where the Jacobian and weight is non
%             %zero
%             cols = ~all(obj.H == 0) & W ~= 0;
%             %Jacobian non zero part
%             H_part = obj.H(:,cols);
%             %Weighted Pinv of the jacobian
%             
%             W_part = diag(W(cols));
%             N = sqrt(W_part);
%             
%             H_part_pinv = pinv(H_part/N);
%             
%             %Kalman Gain Non zero part
%             K_part = K(cols,:);
%             [U,S,V] = svd(K_part);
%             y = V'*obj.dz;
%             z_diff = obj.upper_mes_bound - obj.z_predict.getMesArray()';
%             const_upper = z_diff < 0;
%             Y_upper = U'*H_part_pinv*diag(z_diff);
%             for i=1:size(Y_upper,2)
%                 if(const_upper(i))
%                     %The bound is not sattisfied we need to adjust eigen values 
%                    for j=1:size(S,1)
% 
%                       alpha = Y_upper(j,i)./((S(j,:)*y));
%                       if isinf(alpha) || isnan(alpha) 
%                          alpha = 0; 
%                       end
%                       if alpha > 10
%                           alpha = 10;
%                       end
%                       if alpha < -10
%                           alpha = -10;
%                       end
% 
%                       S(j,:) = S(j,:)*alpha;
%                       alphas(j,i) = alpha;
%                    end       
%                     %Update
%                     K_part = U*S*V';
%                     K(cols,:) = K_part;
%                     x=obj.x_predict+K*obj.dz;                           %Update State Estimate
% 
%                     obj.model_handle.position = x(1:dof);
%                     obj.model_handle.velocity = x(dof+1:end);
%                     obj.model_handle.forwardPosition();
%                     obj.model_handle.forwardVelocity();
%                     z_new = obj.makeMeasure(x).getMesArray()';
%                     
%                     const_upper = (obj.upper_mes_bound-z_new) < 0;
%                     Y_upper = U'*H_part_pinv*diag(obj.upper_mes_bound-z_new);
%                 end
%             end
%             
%             z_diff = obj.z_predict.getMesArray()' - obj.lower_mes_bound;
%             const_lower = z_diff < 0;
%             Y_lower = U'*H_p*diag(z_diff);
%             for i=1:size(Y_lower,2)
%                 if(const_lower(i))
%                     %The bound is not sattisfied we need to adjust eigen values 
%                    for j=1:size(S,1)
% 
%                       alpha = Y_lower(j,i)./((S(j,:)*y));
%                       if isinf(alpha) || isnan(alpha) 
%                          alpha = 1; 
%                       end
%                       if alpha > 1
%                          alpha = 1; 
%                       end
%                       if alpha < -1
%                          alpha = -1; 
%                       end
% 
%                       S(j,:) = S(j,:)*alpha;
%                       alphas(j,i) = alpha;
%                    end   
%                    
%                    %Update 
%                     K = U*S*V';
%                     x=obj.x_predict+K*obj.dz;                           %Update State Estimate
% 
%                     obj.model_handle.position = x(1:dof);
%                     obj.model_handle.velocity = x(dof+1:end);
%                     obj.model_handle.forwardKinematics();
% 
%                      %H = obj.makeH(x);
%                      z_new = obj.makeMeasure(x);
% 
%                      const_lower = (z_new - obj.lower_mes_bound) < 0;
%                      Y_lower = U'*H_p*diag((obj.lower_mes_bound-obj.z_predict));
%                 end
%             end 
        end
        
        function P = makeP(obj)
            
            %Uses Joseph stabalized version of covariance measurement
            %update
            noise = obj.MakeObservationNoise();
            I = eye(size(obj.P_predict));
            P = (I - obj.K*obj.H)*obj.P_predict* ...
                (I - obj.K*obj.H)' + ...
                obj.K*noise(obj.selectedcols,obj.selectedcols)*obj.K';
        end
        
        function Ahat = nearestSPD(obj,A)
            % nearestSPD - the nearest (in Frobenius norm) Symmetric Positive Definite matrix to A
            % usage: Ahat = nearestSPD(A)
            %
            % From Higham: "The nearest symmetric positive semidefinite matrix in the
            % Frobenius norm to an arbitrary real matrix A is shown to be (B + H)/2,
            % where H is the symmetric polar factor of B=(A + A')/2."
            %
            % http://www.sciencedirect.com/science/article/pii/0024379588902236
            %
            % arguments: (input)
            %  A - square matrix, which will be converted to the nearest Symmetric
            %    Positive Definite Matrix.
            %
            % Arguments: (output)
            %  Ahat - The matrix chosen as the nearest SPD matrix to A.

            % test for a square matrix A
            [r,c] = size(A);
            if r ~= c
              error('A must be a square matrix.')
            elseif (r == 1) && (A <= 0)
              % A was scalar and non-positive, so just return eps
              Ahat = eps;
              return
            end

            % symmetrize A into B
            B = (A + A')/2;

            % Compute the symmetric polar factor of B. Call it H.
            % Clearly H is itself SPD.
            [U,Sigma,V] = svd(B);
            H = V*Sigma*V';

            % get Ahat in the above formula
            Ahat = (B+H)/2;

            % ensure symmetry
            Ahat = (Ahat + Ahat')/2;

            % test that Ahat is in fact PD. if it is not so, then tweak it just a bit.
            p = 1;
            k = 0;
            while p ~= 0
              [R,p] = chol(Ahat);
              k = k + 1;
              if p ~= 0
                % Ahat failed the chol test. It must have been just a hair off,
                % due to floating point trash, so it is simplest now just to
                % tweak by adding a tiny multiple of an identity matrix.
                mineig = min(eig(Ahat));
                Ahat = Ahat + (-mineig*k.^2 + eps(mineig))*eye(size(A));
              end
            end
        end
        
    end
end