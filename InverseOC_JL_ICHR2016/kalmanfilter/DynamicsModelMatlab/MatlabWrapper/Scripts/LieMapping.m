classdef LieMapping < handle & matlab.mixin.Copyable & matlab.mixin.Heterogeneous
    %This is EKF implemented with the very basic equations
    
    properties
        %Lie Group Representation
        x_lg;
        dim = [];
        dim_state = [];
    end
    
    properties(SetAccess = protected)
       generators = []; 
    end
    
    properties(Dependent, SetAccess = public)
        %Lie Algebra Representation (hatinv)
        x_la;
        %Lie Algebra Vector Representation (hat)
        x_v;
    end
    
    methods
        % Dependent Setters and getters
        
        function value = get.x_la(obj)
            value = obj.s_log(obj.x_lg);
        end
        
        function set.x_la(obj,value)
            if(-value == value')
                obj.x_lg = obj.s_exp(obj.x_la);
            else
                error('Invalid Algebra');
            end
        end
        
        function set.x_lg(obj,value)
           if(det(value) - 1 < 1e-12)
              obj.x_lg = value;
           else
              error('Invalid Group Element') 
           end
        end
        
        function value = get.x_v(obj)
            value = obj.s_hatinv(obj.x_la);
        end
        
        function set.x_v(obj,value)
            if(size(value,2) == 1)
                obj.x_lg = obj.s_exp(obj.s_hat(value));
            else
                error('Invalid Algebra Vector')
            end
        end
        
        
        % Methods for handling group
        function x_lg = exp(obj)
            %Exponential mapping of an element
            %exp_G : g -> G (LA -> LG)
            x_lg = obj.x_lg;
        end
        
        function x_la = log(obj)
            %Exponential mapping of an element
            %log_G : G -> g (LG -> LA)
            x_la = obj.x_la;
        end
        
        function x_v = hatinv(obj)
            %[.]¡_G : g -> R^p
            %x_la expected as member of Lie Algebra joint to G
            x_v = obj.x_v;
        end
        
       
        function adj_x = adjointSmall(obj,x_v)
            %ad_G : g -> L(g)
            %x_v expected as a vector
            adj_x = zeros(LieMapping.dim,1);
        end
        
        function Adj_x = adjointLarge(obj,x_lg)
            %Ad_G : G -> GL(g)
            %x_lg expected as member of Lie Group G
            Adj_x = zeros(LieMapping.dim,1);
        end 
    end
    
    methods (Static)
        %Declaration of static LieMapping functions 
        function x_la = s_hat(x_v)
            %Representation adaptation of an element
            %[.]^_G : R^p -> g
            %x_v expected as a vector
            x_la = [];
        end
        
        function x_lg = s_exp(x_la)
            %Exponential mapping of an element
            %exp_G : g -> G (LA -> LG)
            %x_la expected as member of Lie Algebra joint to G
            x_lg = [];
        end
        
        function x_la = s_log(x_lg)
            %Exponential mapping of an element
            %log_G : G -> g (LG -> LA)
            %x_lg expected as member of Lie Group G
            x_la = [];
        end
        
        function x_v = s_hatinv(x_la)
            x_v = [];
        end
        
    end
    
end

