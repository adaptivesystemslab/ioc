function [ceq, ceq_x, ceq_dx, ceq_ddx] = calc_direct_const(c_const, feature_use, param)
    if isfield(param, 'const_x') && ~isempty(param.const_x)
        const_y = sym(param.const_y);
        ceq_x = calcConstraint(feature_use.q, param.const_x, const_y, c_const(1));
    else
        ceq_x = [];
    end

    if isfield(param, 'const_dx') && ~isempty(param.const_dx)
        const_dy = sym(param.const_dy);
        ceq_dx = calcConstraint(feature_use.dq, param.const_dx, const_dy, c_const(2));
    else
        ceq_dx = [];
    end

    if isfield(param, 'const_ddx') && ~isempty(param.const_ddx)
        const_ddy = sym(param.const_ddy);
        ceq_ddx = calcConstraint(feature_use.ddq, param.const_ddx, const_ddy, c_const(3));
    else
        ceq_ddx = [];
    end

    ceq = [ceq_x; ceq_dx; ceq_ddx];
end

function ceq = calcConstraint(featureVal, const_x, const_y, c_const)
    ceqTemp = sym(zeros(size(const_y)));
    
    % look at each timestamp
    for jj = 1:length(const_x) 
        ceqTemp(:, jj) = featureVal(:, const_x(jj)) - const_y(:, jj);
    end
    
    ceq = c_const*mat2vec(ceqTemp);
    
    % look at each dof
%     for jj = 1:size(featureVal, 1)
%         temp = featureVal(jj, const_x) - const_y(jj, :);
%         ceqTemp(jj, :) = temp;
% %         ceqTemp(jj, :) = temp / max(max(abs(temp)));
%     end
%     ceq = c_const*ceqTemp(:);
end