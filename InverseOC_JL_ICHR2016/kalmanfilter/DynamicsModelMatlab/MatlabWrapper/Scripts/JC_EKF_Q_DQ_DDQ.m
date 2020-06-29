classdef JC_EKF_Q_DQ_DDQ < EKF_Q_DQ_DDQ
    %This is EKF that uses a Kinematic Model with sensors and can have
    %constraints on state and measurement using Gain Projection
    
    properties     
       upper_state_bound = [];
       lower_state_bound = [];
       upper_mes_bound = [];
       lower_mes_bound = [];
    end
    
    methods
        function obj = JC_EKF_Q_DQ_DDQ(model)
           obj = obj@EKF_Q_DQ_DDQ(model);
           
           %Init Bounds
           obj.upper_state_bound = ones(obj.sizeX,1).*inf;
           obj.lower_state_bound = ones(obj.sizeX,1).*-inf;
           
           obj.upper_mes_bound = ones(obj.sizeZ,1).*inf;
           obj.lower_mes_bound = ones(obj.sizeZ,1).*-inf;
        end
        
        
%         function run_iteration(obj,u,z,matches)
%             obj.run_iteration@EKF_Resize(u,z,matches);
%             obj.makeH(obj.state);
%             
%             J = obj.model_handle.sensors(2).baseJacobian;
%             Jd = obj.model_handle.sensors(2).baseJacobianDerivative;
%             D = [zeros(6,size(J,2)) J zeros(6,size(J,2))];
%             obj.state=obj.state - D'/(D*D'+1e-12*eye(size(D,1)))*D*obj.state;
%         end
        
    end
    
    methods(Access=protected)
        
        function K = makeK(obj)
            %Call regular make K
            K = obj.makeK@EKF();
            %Predict
            x = obj.x_predict+K*obj.dz;
            obj.makeH(x);
            J = obj.model_handle.sensors(2).baseJacobian;
            Jd = obj.model_handle.sensors(2).baseJacobianDerivative;
            D = [zeros(6,size(J,2)) J zeros(6,size(J,2));...
                zeros(6,size(J,2)) Jd J];
            c = cond(D*D'+1e-8*eye(size(D,1)))
            if c < 300
                K = K - D'/(D*D'+1e-8*eye(size(D,1)))*(D*x)/(obj.dz'/obj.inov_covariance*obj.dz)*obj.dz'/obj.inov_covariance;
            end
        end
        
        function P = makeP(obj)
          P = obj.P_predict-obj.K*obj.H*obj.P_predict;
          P = obj.nearestSPD(P);
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