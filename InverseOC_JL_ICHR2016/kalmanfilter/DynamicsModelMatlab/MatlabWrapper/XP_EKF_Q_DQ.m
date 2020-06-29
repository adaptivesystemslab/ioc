classdef XP_EKF_Q_DQ < EKF_Q_DQ
    %This is EKF that uses a Kinematic Model with sensors and can have
    %constraints on state and measurement using SVD
    
    properties
        upper_state_bound = [];
        lower_state_bound = [];
        upper_mes_bound = [];
        lower_mes_bound = [];
    end
    
    methods
        function obj = XP_EKF_Q_DQ(model)
            obj = obj@EKF_Q_DQ(model);
            
            %Init Bounds
            obj.upper_state_bound = ones(obj.sizeX,1).*inf;
            obj.lower_state_bound = ones(obj.sizeX,1).*-inf;
            
            obj.upper_mes_bound = ones(obj.sizeZ,1).*inf;
            obj.lower_mes_bound = ones(obj.sizeZ,1).*-inf;
        end
    end
    
    methods(Access=protected)
        
        function P = makeP(obj)
            
            %Uses Joseph stabalized version of covariance measurement
            %update
            noise = obj.MakeObservationNoise();
            I = eye(size(obj.P_predict));
            P = (I - obj.K*obj.H)*obj.P_predict* ...
                (I - obj.K*obj.H)' + ...
                obj.K*noise(obj.selectedcols,obj.selectedcols)*obj.K';
        end
        
        function x = makeState(obj)
            %Project state onto constraints
            x_orig = obj.makeState@EKF_Q_DQ();
            
            const_upper = x_orig > obj.upper_state_bound;
            %This is q dqformulation of EKF, if q constraint is not
            %sattisfied impose constraint on dq to be 0
            dof = numel(obj.model_handle.joints);
            constraints = x_orig;
            constraints(const_upper) = obj.upper_state_bound(const_upper);
            I = find(const_upper(1:dof));
            if ~isempty(I)
                %const_upper(I+dof) = x(I+dof) > 0;
                constraints(I+dof) = 0;
            end
            
            const_lower = x_orig < obj.lower_state_bound;
            constraints(const_lower) = obj.lower_state_bound(const_lower);
            I = find(const_lower(1:dof));
            if ~isempty(I)
                %const_lower(I+dof) = x(I+dof) > 0;
                constraints(I+dof) = 0;
            end
            
            d = constraints;
            %Project Gain
            K = K - D'/(D*D')*(D*x-d)/(obj.dz'/obj.inov_covariance*obj.dz)*obj.dz'/obj.inov_covariance;
            
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