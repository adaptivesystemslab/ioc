classdef Rn < LieMapping
    %This is EKF implemented with the very basic equations
    
    properties(Constant)
        
    end
    
    methods
        % Methods for handling group
        
        function obj = Rn(value)
           %Rn constructor 
           if size(value,2) == 1
               %We have a vector input
               obj.x_v = value;
           elseif size(value,1) == size(value,2) &...
                   diag(value) == 0
               %We have a lie algebra input
               obj.x_la = value;
           elseif size(value,1) == size(value,2)
               %We have a lie group input
               obj.x_lg = value;
           else
               error('Unknown RN format');
           end
           obj.dim = numel(obj.x_v);
           obj.dim_state = size(obj.x_lg,1);
           
           %Create generators
           base_gen = zeros(obj.dim_state);
           for i=1:obj.dim
              obj.generators(:,:,i) = base_gen;
              obj.generators(i,end,i) = 1;
           end
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
            adj_x = zeros(numel(obj.x_v));
        end
        
        function Adj_x = adjointLarge(obj)
            %Ad_G : G -> GL(g)
            %x_lg expected as member of Lie Group G
            Adj_x = obj.s_adjointLarge(obj.x_lg);
        end
        
    end
    
    methods(Static)
        function x_lg = s_exp(x_la)
            %Exponential mapping of an element
            %exp_G : g -> G (LA -> LG)
            %x_la expected as member of Lie Algebra joint to G
            x_lg = eye(size(x_la,1)) + x_la;
        end
        
        function x_la = s_log(x_lg)
            %Exponential mapping of an element
            %log_G : G -> g (LG -> LA)
            %x_lg expected as member of Lie Group G
            x_la = zeros(size(x_lg));
            x_la(1:end-1,end) = x_lg(1:end-1,end);
        end
        
        function x_la = s_hat(x_v)
            %Representation adaptation of an element
            %[.]^_G : R^p -> g
            %x_v expected as a vector
            x_la = zeros(length(x_v)+1);
            x_la(1:length(x_v),end) = x_v(1:end);
        end
        
        function x_v = s_hatinv(x_la)
            %[.]¡_G : g -> R^p
            %x_la expected as member of Lie Algebra joint to G
            x_v = x_la(1:end-1,end);
        end
        
        function adj_x = s_adjointSmall(x_v)
            %ad_G : g -> L(g)
            %x_v expected as a vector
            adj_x = zeros(numel(x_v));
        end
        
        function Adj_x = s_adjointLarge(x_lg)
            %Ad_G : G -> GL(g)
            %x_lg expected as member of Lie Group G
            Adj_x = eye(size(x_lg,1)-1);
        end

    end
    
end