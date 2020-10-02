classdef SVD_QP_EKF_Q_DQ < EKF_Q_DQ
    %This is EKF that uses a Kinematic Model with sensors and can have
    %constraints on state and measurement using SVD
    
    properties     
       upper_state_bound = [];
       lower_state_bound = [];
       upper_mes_bound = [];
       lower_mes_bound = [];
    end
    
    methods
        function obj = SVD_QP_EKF_Q_DQ(model)
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
            
            
            if isempty(find(obj.z_predict.getMesArray' + (obj.H*K*obj.dz)  > obj.upper_mes_bound,1))
               return; 
            end
            
            W = [1 1 1 0 0 0];
            
            %Pull out only columns where the Jacobian and weight is non
            %zero
            cols = ~all(obj.H == 0) & W ~= 0;
            %Jacobian non zero part
            H_part = obj.H(:,cols);
            W_part = diag(W(cols));
            K_part = K(cols,:);
            K_vec = vec(K_part);
            e_x = x(cols) - obj.x_predict(cols);
            
            %Set up the quad prog
            dz_I_kron = kron(obj.dz',eye(size(K_part,1)));
            H = dz_I_kron'*W_part*dz_I_kron;
            H = H+eye(size(H))*1e-7;
            f = kron(obj.dz',e_x'*W_part)';
            
            valid_bnd_indxs = obj.upper_mes_bound < inf;
            z_pred_array = obj.z_predict.getMesArray();
            z_pred_array = z_pred_array(valid_bnd_indxs);
            b = obj.upper_mes_bound(valid_bnd_indxs) - z_pred_array;
            A = kron(obj.dz',H_part);
            A = A(valid_bnd_indxs,:);
            
            K_part_prime_vec = quadprog(H,f,A,b,[],[],[],[],vec(K_part));
            K_part_prime  = vec_inv(K_part_prime_vec,size(K_part,1));
            
            K(cols,:) = K_part_prime;
            
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