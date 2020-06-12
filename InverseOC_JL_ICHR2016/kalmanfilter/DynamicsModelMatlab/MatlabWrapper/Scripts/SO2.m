classdef SO2 < LieMapping
    %SO2 Element
    
    properties
        
    end
    
    methods
        % Methods for handling group
        
        function obj = SO2(value)
            %Constructor of SO2 element
            %The input can be exp or log
            
            if numel(value) == 1
                %Input is lie algebra vector representation
                obj.x_v = value;
            elseif (-value == value')
                %Input is Lie Algebra representation
                obj.x_la = value;
            elseif (size(value,1) == 2 && size(value,2) == 2 && det(value)-1 < 1e-10)
                %Input is a Lie Group Element
                obj.x_lg = value;
            else
                error('Incorrect SO2 Constructor Input');
            end
            obj.dim = 1;
            obj.dim_state = 2; % means stored in 2x2 matrices
        end
        
        function adj_x = adjointSmall(obj)
            %ad_G : g -> L(g)
            %x_v expected as a vector
            adj_x = zeros(1);
        end
        
        function Adj_x = adjointLarge(obj)
            %Ad_G : G -> GL(g)
            %x_lg expected as member of Lie Group G
            Adj_x = eye(1);
        end
        
    end
    
    methods(Static)
        function x_la = s_hat(x_v)
            %Representation adaptation of an element
            %[.]^_G : R^p -> g
            %x_v expected as a vector
            x_la = [...
                0 -x_v ;...
                x_v 0];
        end
        
        function x_lg = s_exp(x_la)
            %Exponential mapping of an element
            %exp_G : g -> G (LA -> LG)
            %x_la expected as member of Lie Algebra joint to G
            fi = x_la(2,1);
            x_lg = cos(fi) * eye(2) + sin(fi) * [0 -1; 1 0];
        end
        
        function x_la = s_log(x_lg)
            %Exponential mapping of an element
            %log_G : G -> g (LG -> LA)
            %x_lg expected as member of Lie Group G
            fi = atan2(x_lg(2,1), x_lg(1,1));
            x_la = fi * [0 -1; 1 0];
        end
        
        function x_v = s_hatinv(x_la)
            x_v = x_la(2,1);
        end
        
        function adj_x = s_adjointSmall(x_v)
            %ad_G : g -> L(g)
            %x_v expected as a vector
            adj_x = zeros(1);
        end
        
        function Adj_x = s_adjointLarge(x_lg)
            %Ad_G : G -> GL(g)
            %x_lg expected as member of Lie Group G
            Adj_x = eye(1);
        end
    end
    
end