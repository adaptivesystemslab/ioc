classdef SE3 < LieMapping
    %This is EKF implemented with the very basic equations
    
    properties
        
    end
    
    methods
        % Methods for handling group
        
        function obj = SE3(value)
            %Constructor of SE3 object
            %Input can be:
            %SO2 to upcast into SE3 assuming rot about Z and no translation
            %SO3 to upcast into SE3 assuming no translation
            %Lie Group element (Transformation Matrix)
            %Lie Algebra element
            %Lie Algebra vector element
            if nargin ~= 0
                %For Upcasting from SO2
                %Check if we were passed SO2 or SO3 objects
                if(isa(value,'LieMapping'))
                    [m,n] = size(value);
                    for i=1:m
                        for j=1:n
                            obj(i,j) = SE3();
                            cur_obj = value(i,j);
                            if isa(value(i,j),'SO2')
                                obj(i,j).x_lg = [cur_obj.x_lg [0 0; 0 0]; 0 0 1 0; 0 0 0 1];
                                %Only single generator in SE3 format
                                obj(i,j).generators = [[0 -1 0 0; 1 0  0 0;  0 0 0 0]; zeros(1,3) 0];
                                obj(i,j).dim = 6;
                                obj(i,j).dim_state = 4; % means stored in 4x4 matrices
                                %For Upcasting from SO3
                            elseif isa(value(i,j),'SO3')
                                obj(i,j).x_lg = [cur_obj.x_lg [0;0;0]; 0 0 0 1];
                                %Rotation Generators when upcasting from SO3
                                obj(i,j).generators = zeros(4,4,3);
                                obj(i,j).generators(:,:,1) = [[0 0  0 0; 0 0 -1 0;  0 1 0 0]; zeros(1,3) 0];
                                obj(i,j).generators(:,:,2) = [[0 0  1 0; 0 0  0 0; -1 0 0 0]; zeros(1,3) 0];
                                obj(i,j).generators(:,:,3) = [[0 -1 0 0; 1 0  0 0;  0 0 0 0]; zeros(1,3) 0];
                                obj(i,j).dim = 6;
                                obj(i,j).dim_state = 4; % means stored in 4x4 matrices
                            elseif isa(value(i,j),'SE3')
                                obj(i,j).x_lg = cur_obj.x_lg;
                                %Rotation Generators when upcasting from SO3
                                obj(:,:).generators = zeros(4,4,6);
                                obj(:,:).generators(:,:,1) = [[0 0  0 0; 0 0 -1 0;  0 1 0 0]; zeros(1,3) 0];
                                obj(:,:).generators(:,:,2) = [[0 0  1 0; 0 0  0 0; -1 0 0 0]; zeros(1,3) 0];
                                obj(:,:).generators(:,:,3) = [[0 -1 0 0; 1 0  0 0;  0 0 0 0]; zeros(1,3) 0]; % around Z????
                                obj(:,:).generators(:,:,4) = [[0 0  0 1; 0 0 0 0;  0 0 0 0]; zeros(1,3) 0];
                                obj(:,:).generators(:,:,5) = [[0 0  0 0; 0 0 0 1;  0 0 0 0]; zeros(1,3) 0];
                                obj(:,:).generators(:,:,6) = [[0 0  0 0; 0 0 0 0;  0 0 0 1]; zeros(1,3) 0];
                                obj(:,:).dim = 6;
                                obj(:,:).dim_state = 4; % means stored in 4x4 matrices
                            end
                        end
                    end
                elseif size(value,1) == 6 && size(value,2) == 1
                    %Vector input, must be x_v
                    obj.x_v = value;
                elseif -value == value'
                    %input is algebra be Algebra element
                    obj.x_la = value;
                elseif det(value(1:3,1:3))-1 < 1e-10 && value(4,4) == 1
                    obj.x_lg = value;
                else
                    error('Incorrect SE3 Constructor Input');
                end
                if(isempty(obj(1).generators))
                    obj(:,:).generators = zeros(4,4,6);
                    obj(:,:).generators(:,:,1) = [[0 0  0 0; 0 0 -1 0;  0 1 0 0]; zeros(1,3) 0];
                    obj(:,:).generators(:,:,2) = [[0 0  1 0; 0 0  0 0; -1 0 0 0]; zeros(1,3) 0];
                    obj(:,:).generators(:,:,3) = [[0 -1 0 0; 1 0  0 0;  0 0 0 0]; zeros(1,3) 0]; % around Z????
                    obj(:,:).generators(:,:,4) = [[0 0  0 1; 0 0 0 0;  0 0 0 0]; zeros(1,3) 0];
                    obj(:,:).generators(:,:,5) = [[0 0  0 0; 0 0 0 1;  0 0 0 0]; zeros(1,3) 0];
                    obj(:,:).generators(:,:,6) = [[0 0  0 0; 0 0 0 0;  0 0 0 1]; zeros(1,3) 0];
                    obj(:,:).dim = 6;
                    obj(:,:).dim_state = 4; % means stored in 4x4 matrices
                end
            end
        end
        
        function adj_x = adjointSmall(obj)
            %ad_G : g -> L(g)
            %x_v expected as a vector
            x_v = obj.x_v;
            adj_x = [tmp.hat(x_v(1:3)) SO3.s_hat(x_v(4:6)); zeros(3) SO3.s_hat(x_v(1:3))];
        end
        
        function Adj_x = adjointLarge(obj)
            %Ad_G : G -> GL(g)
            %x_lg expected as member of Lie Group G
            x_lg = obj.x_lg;
            Adj_x = [x_lg(1:3,1:3) SO3.s_hat(x_lg(1:3,4))*x_lg(1:3,1:3); zeros(3) x_lg(1:3,1:3)];
        end
    end
    
    methods(Static)
        function x_lg = s_exp(x_la)
            %Exponential mapping of an element
            %exp_G : g -> G (LA -> LG)
            %x_la expected as member of Lie Algebra joint to G
            
            omega = [x_la(3,2);x_la(1,3);x_la(2,1)];
            theta = norm(omega);
            epxp_omega = eye(3);
            V = eye(3);
            if theta ~= 0
                epxp_omega = epxp_omega + sin(theta)/theta*skew(omega) + ...
                    (1-cos(theta))/theta^2*skew(omega)^2;
                V = V + (1-cos(theta))/theta^2*skew(omega) + ...
                    (theta-sin(theta))/theta^3*skew(omega)^2;
            end
            x_lg = [epxp_omega V*x_la(1:3,4); zeros(1,3) 1];
            
        end
        
        function x_la = s_log(x_lg)
            %Exponential mapping of an element
            %log_G : G -> g (LG -> LA)
            %x_lg expected as member of Lie Group G
            
            R = x_lg(1:3,1:3);
            theta = acos((trace(R)-1)/2);
            ln_R = theta/(2*sin(theta))*(R-R');
            if theta ~= 0
                V = eye(3) + (1-cos(theta))/theta^2*ln_R + (theta-sin(theta))/theta^3*ln_R^2;
            else
                ln_R = zeros(3);
                V = eye(3);
            end
            x_la = [ln_R V\x_lg(1:3,4); zeros(1,3) 0];
        end
        
        function x_la = s_hat(x_v)
            %Representation adaptation of an element
            %[.]^_G : R^p -> g
            %x_v expected as a vector
            x_rot = SO3(x_v(1:3)).x_la;
            x_la = [x_rot x_v(4:6); zeros(1,3) 0];
        end
        
        function x_v = s_hatinv(x_la)
            %[.]¡_G : g -> R^p
            %x_la expected as member of Lie Algebra joint to G
            x_v = [x_la(3,2) x_la(1,3) x_la(2,1) x_la(1:3,4)']';
        end
        
        function adj_x = s_adjointSmall(x_v)
            %ad_G : g -> L(g)
            %x_v expected as a vector
            adj_x = [SO3.s_hat(x_v(1:3)) SO3.s_hat(x_v(4:6)); zeros(3) SO3.s_hat(x_v(1:3))];
        end
        
        function Adj_x = s_adjointLarge(x_lg)
            %Ad_G : G -> GL(g)
            %x_lg expected as member of Lie Group G
            Adj_x = [x_lg(1:3,1:3) SO3.s_hat(x_lg(1:3,4))*x_lg(1:3,1:3); zeros(3) x_lg(1:3,1:3)];
        end
    end
    
    methods(Static)
        function invT = fastInverse(T)
            invT = T';
            invT(4,1:3) = 0;
            invT(1:3,4) = -invT(1:3,1:3)*T(1:3,4);
        end
    end
    
end