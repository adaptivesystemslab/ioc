classdef SO3 < LieMapping
    %This is EKF implemented with the very basic equations
    
    properties
        
    end
    
    methods
        % Methods for handling group
        
        function obj = SO3(value)
            %Constructor of SO3 object
            %Input can be:
            %SO2 to upcast into SO3 assuming rotation about Z
            %SE3 to downcast into SO3 assuming zero translation
            %Lie Group element (Rotation Matrix)
            %Lie Algebra element
            %Lie Algebra vector element
            
            %For Upcasting SO2
            if isa(value,'SO2')
                value = [value.x_lg [0; 0]; 0 0 1];
            %For downcasting into SO3 from SE3 if translational component
            %is 0
            elseif isa(value,'SE3') && value.x_lg(1:3,4) == 0
                value = value.x_lg(1:3,1:3);
            end
                
            if size(value,1) == 3 && size(value,2) == 1
               %Vector input, must be x_v
               obj.x_v = value;
            elseif -value == value' 
                %input is algebra be Algebra element
                obj.x_la = value;
            elseif size(value,1) == 3 && size(value,2) == 3 && ...
                    det(value)-1 < 1e-10
                obj.x_lg = value;
            else
                error('Incorrect SO3 Constructor Input');
            end
            
            obj.dim = 3;
            obj.dim_state = 3; % means stored in 3x3 matrices
            
            obj.generators = zeros(3,3,3);
            obj.generators(:,:,1) = [0 0  0 ; 0 0 -1 ;  0 1 0 ];
            obj.generators(:,:,2) = [0 0  1 ; 0 0  0 ; -1 0 0 ];
            obj.generators(:,:,3) = [0 -1 0 ; 1 0  0 ;  0 0 0 ];
        end
        
        function x_lg = exp(obj)
            %Exponential mapping of an element
            %exp_G : g -> G (LA -> LG)
            x_lg = obj.x_lg;
        end
        
        function x_la = log(obj)
            %Exponential mapping of an element
            %log_G : G -> g (LG -> LA)
            %x_lg expected as member of Lie Group G
            x_la = obj.x_la;
        end
        
        function x_la = hat(obj)
            %Representation adaptation of an element
            %[.]^_G : R^p -> g
            %x_v expected as a vector
            x_la = obj.x_la;
        end
        
        function x_v = hatinv(obj)
            %[.]¡_G : g -> R^p
            x_v = obj.x_v;
        end
        
        function adj_x = adjointSmall(obj)
            %ad_G : g -> L(g)
            adj_x = obj.x_la;
        end
        
        function Adj_x = adjointLarge(obj)
            %Ad_G : G -> GL(g)
            %x_lg expected as member of Lie Group G
            Adj_x = obj.x_lg;
        end
        
    end
    
    methods(Static)
       function x_lg = s_exp(x_la)
            %Exponential mapping of an element
            %exp_G : g -> G (LA -> LG)
            %x_la expected as member of Lie Algebra joint to G
            phi = [x_la(3,2);x_la(1,3);x_la(2,1)];
            nphi = norm(phi);
            min = 10^-10;
            % NaNs appear if norm(phi)==0
            if norm(phi)==0
                x_la(2,1) = min; x_la(3,2) = min; x_la(1,3) = min;
                x_la(1,2) = -min; x_la(2,3) = -min; x_la(3,1) = -min;
                phi = [x_la(3,2);x_la(1,3);x_la(2,1)];
                nphi = mod(norm(phi),pi);
            end
            % prema Barfoot-u
            x_lg = cos(nphi)*eye(3) + sin(nphi)/nphi * x_la + (1-cos(nphi))/nphi^2 * phi * phi'; %% CHECK however!!!!
            %x_lg = eye(3) + sin(nphi)/nphi * x_la + (1-cos(nphi))/nphi^2 * (x_la * x_la);
        end
        
        function x_la = s_log(x_lg)
            %Exponential mapping of an element
            %log_G : G -> g (LG -> LA)
            %x_lg expected as member of Lie Group G
            trR = trace(x_lg);
            w = acos((trR - 1)/2);
            if w==0
                x_la = zeros(3);
            else
                x_la = w/(2*sin(w)) * (x_lg-x_lg');
            end
        end
        
        function x_la = s_hat(x_v)
            %Representation adaptation of an element
            %[.]^_G : R^p -> g
            %x_v expected as a vector
            x_la = [...
                0        -x_v(3)  x_v(2) ;...
                x_v(3)   0        -x_v(1) ;...
                -x_v(2)  x_v(1)   0];
        end
        
        function x_v = s_hatinv(x_la)
            %[.]¡_G : g -> R^p
            %x_la expected as member of Lie Algebra joint to G
            x_v = [x_la(3,2) x_la(1,3) x_la(2,1)]';
        end
        
        function adj_x = s_adjointSmall(x_v)
            %ad_G : g -> L(g)
            %x_v expected as a vector
            adj_x = SO3.s_hat(x_v(1:3));
        end
        
        function Adj_x = s_adjointLarge(x_lg)
            %Ad_G : G -> GL(g)
            %x_lg expected as member of Lie Group G
            Adj_x = x_lg;
        end 
    end
end