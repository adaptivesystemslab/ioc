function [ceq, ceq_x, ceq_dx, ceq_ddx, ceq_dyn] = calc_direct_const(c_const, feature_use, param)
    if isfield(param, 'const_x') && ~isempty(param.const_x)
        ceq_x = calcConstraint(feature_use.q, param.const_x, param.const_y, c_const(1));
    end

    if isfield(param, 'const_dx') && ~isempty(param.const_dx)
        ceq_dx = calcConstraint(feature_use.dq, param.const_dx, param.const_dy, c_const(2));
    else
        ceq_dx = [];
    end

    if isfield(param, 'const_ddx') && ~isempty(param.const_ddx)
        ceq_ddx = calcConstraint(feature_use.ddq, param.const_ddx, param.const_ddy, c_const(3));
    else
        ceq_ddx = [];
    end
    
%     gapT = sort([param.const_x-1 param.const_x param.const_x+1]);
%     ceq_dyn = calcRelations(feature_use, gapT, c_const(4:6));

%     ceq = [ceq_x; ceq_dx; ceq_ddx; ceq_dyn];
    ceq = [ceq_x; ceq_dx; ceq_ddx];
end

% function [ceq, ceq_x, ceq_dx, ceq_ddx, ceq_dyn] = calc_direct_const(c_const, feature_use, param)
%     if isfield(param, 'const_x') && ~isempty(param.const_x)
%         ceq_x = calcConstraint(feature_use.q, param.const_x, param.const_y, c_const(1));
%     else
%         ceq_x = [];
%     end
% 
%     if isfield(param, 'const_dx') && ~isempty(param.const_dx)
%         ceq_dx = calcConstraint(feature_use.dq, param.const_dx, param.const_dy, c_const(2));
%     else
%         ceq_dx = [];
%     end
% 
%     if isfield(param, 'const_ddx') && ~isempty(param.const_ddx)
%         ceq_ddx = calcConstraint(feature_use.ddq, param.const_ddx, param.const_ddy, c_const(3));
%     else
%         ceq_ddx = [];
%     end
%     
%     gapT = sort([param.const_x-1 param.const_x param.const_x+1]);
%     ceq_dyn = calcRelations(feature_use, gapT, c_const(4:6));
% 
%     ceq = [ceq_x; ceq_dx; ceq_ddx; ceq_dyn];
% end

function ceq = calcConstraint(featureVal, const_x, const_y, c_const)
    ceqTemp = zeros(size(const_y));
    
    % look at each timestamp
    for jj = 1:length(const_x) 
        ceqTemp(:, jj) = featureVal(:, const_x(jj)) - const_y(:, jj);
    end
    
    ceq = c_const*mat2vec(ceqTemp);
end

function ceq = calcRelations(featureVal, t_skip, c_const)
%     ceqTemp = zeros(size(const_y));
    
    % look at each timestamp
    ceq_ddq = zeros(size(featureVal.ddq, 1), size(featureVal.ddq_tp1, 2) - size(t_skip, 2));
    counter = 0;
    for jj = 2:size(featureVal.ddq_tp1, 2)
        if sum(jj == t_skip) == 0
           counter = counter + 1; 
           ceq_ddq(:, counter) = featureVal.ddq_tp1(:, jj-1) - featureVal.ddq(:, jj);
        end
        
%         ceq_dq(:, jj)  = featureVal.dq_tp1(:, jj)  - featureVal.dq(:, jj+1);
%         ceq_q(:, jj)   = featureVal.q_tp1(:, jj)   - featureVal.q(:, jj+1);
    end
    
%     ceqTemp = [c_const(1)*ceq_ddq; c_const(2)*ceq_dq; c_const(3)*ceq_q];
    ceqTemp = [c_const(1)*ceq_ddq;];
    ceq = mat2vec(ceqTemp);
end