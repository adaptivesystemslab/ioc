function [ceq_outlin, ceq_x, ceq_dx, ceq_ddx] = calc_direct_const(c_const, feature_use, param)
    ceq_x = zeros(size(feature_use.q));
    ceq_dx = zeros(size(feature_use.dq));
    ceq_ddx = zeros(size(feature_use.ddq));
    lenT = size(feature_use.q, 2);
    
    if isfield(param, 'const_x') && ~isempty(param.const_x)
        ceq_x = calcConstraint(ceq_x, feature_use.q, param.const_x, param.const_y, c_const(1));
        x_index = param.const_x;
    end

    if isfield(param, 'const_dx') && ~isempty(param.const_dx)
        ceq_dx = calcConstraint(ceq_dx, feature_use.dq, param.const_dx, param.const_dy, c_const(2));
        dx_index = param.const_x;
    end

    if isfield(param, 'const_ddx') && ~isempty(param.const_ddx)
        ceq_ddx = calcConstraint(ceq_ddx, feature_use.ddq, param.const_ddx, param.const_ddy, c_const(3));
        gapT = sort([param.const_x]);
%         gapT = [];
    end
    
    deltaT = param.dynamicSampleRate;
%     applicationZone = [deltaT:deltaT:(lenT/2-deltaT) (lenT/2+deltaT):(lenT-deltaT)];
    applicationZone = 1:deltaT:lenT;
    if deltaT < lenT && deltaT > 0
        % if the sampling rate is appropriate
        
        ceq_ddx = calcRelations(ceq_ddx, feature_use, gapT, c_const(4:6), applicationZone);
        ddx_index = unique([param.const_x applicationZone]);
    else
        % otherwise, assume dynamic update is off
        ddx_index = param.const_x;
    end
    
    % remove excess zeros
    ceq = [ceq_x(:, x_index) ceq_dx(:, dx_index) ceq_ddx(:, ddx_index)];
    ceq_outlin = mat2vec(ceq);
end

function ceqTemp = calcConstraint(ceqTemp, featureVal, const_x, const_y, c_const)
    % look at each timestamp
    for jj = 1:length(const_x) 
        ceqTemp(:, const_x(jj)) = c_const*(featureVal(:, const_x(jj)) - const_y(:, jj));
    end
end

function ceqTemp = calcRelations(ceqTemp, featureVal, t_skip, c_const, applicationZone)
%     ceqTemp = zeros(size(const_y));
    
    % look at each timestamp
%     ceq_ddq = zeros(size(featureVal.ddq, 1), size(featureVal.ddq_tp1, 2) - size(t_skip, 2));
%     counter = 0;
    for jj = applicationZone
        if sum(jj == t_skip) == 0
           ceqTemp(:, jj) = c_const(1)*(featureVal.ddq_tp1(:, jj-1) - featureVal.ddq(:, jj));
        else
            llala = 0;
        end
        
%         ceq_dq(:, jj)  = featureVal.dq_tp1(:, jj)  - featureVal.dq(:, jj+1);
%         ceq_q(:, jj)   = featureVal.q_tp1(:, jj)   - featureVal.q(:, jj+1);
    end
    
% %     ceqTemp = [c_const(1)*ceq_ddq; c_const(2)*ceq_dq; c_const(3)*ceq_q];
%     ceqTemp = [c_const(1)*ceq_ddq;];
%     ceq = mat2vec(ceqTemp);
end