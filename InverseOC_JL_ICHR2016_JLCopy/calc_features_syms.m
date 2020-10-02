function feature_out = calc_features_syms(splineFit, feature_win, param)
    % given the spline variables, produce all the individual variables
        
    % calculate the joint angles
    if isstruct(splineFit) 
        if isfield(splineFit, 'syms') && splineFit.syms
            q = sym('q', [param.dof_count param.NbSample]);
            dq = sym('dq', [param.dof_count param.NbSample]);
            ddq = sym('ddq', [param.dof_count param.NbSample]);
            dddq = sym('dddq', [param.dof_count param.NbSample]);
            
            [feature_out.q_syms, feature_out.dq_syms, feature_out.ddq_syms, feature_out.dddq_syms, feature_out.t_syms_range] ...
                = calc_q_syms(splineFit.q, param);

            fctMode = 'symbolic';
        else
            q = splineFit.q;
            dq = splineFit.dq;
            ddq = splineFit.ddq;
            dddq = splineFit.dddq;
            
            fctMode = 'straight';
        end
    else
        % spline fit is actually a q array
        q = splineFit;
        
        dq = calcDeriv(q, param.dt_spline);
        ddq = calcDeriv(dq, param.dt_spline);
        dddq = calcDeriv(ddq, param.dt_spline);
        
        fctMode = 'straight';
    end
    
    % if the 'param' internal time length and the 'splineFit' internal time
    % length is not the same, update 'param' one so that the dynamic
    % calculations will match up properly
    if length(param.t_spline) ~= size(q, 2)
        param = update_length(param, size(q, 2));
    end
    
    % now need to maniuplate q matrix to have 7DOFs
    lengthQ = size(q, 2);
    q_7dof   = [zeros(param.dofsFromFull(1)-1, lengthQ);   q;   zeros(param.dofsFull(end)-param.dofsFromFull(end), lengthQ)];
    dq_7dof  = [zeros(param.dofsFromFull(1)-1, lengthQ);   dq;  zeros(param.dofsFull(end)-param.dofsFromFull(end), lengthQ)];
    ddq_7dof = [zeros(param.dofsFromFull(1)-1, lengthQ);   ddq; zeros(param.dofsFull(end)-param.dofsFromFull(end), lengthQ)];
    z_7dof = zeros(size(q_7dof));
    
%     % if q is shorter than param, then pad q_(t=n)s at the end
%     if size(q, 2) < param.NbSample + 6 % +4 for the dddq
%         q = q(:, [1 1 1 1:end end end end]);
%     end
%     
    % then calculate the dq-dddq
%     dq = calcDeriv(q, param.dt_spline);
%     ddq = calcDeriv(dq, param.dt_spline);
%     dddq = calcDeriv(ddq, param.dt_spline);
    
%     % now crop the q-dddq so it's back down to proper length
%     newInd = 4:(size(q, 2)-3);
%     q = q(:, newInd);
%     dq = dq(:, newInd);
%     ddq = ddq(:, newInd);
%     dddq = dddq(:, newInd);
    
    % now calculate the forward kinematics and obtain EF pos/vel/acc
    init_mat = repmat(eye(4),[1 1 param.NbSample]);
    Tr = FKM_SQUAT_7DOF(init_mat, q_7dof', param.L'); % end eff 
    J = SQUAT_7DOF_Ext_wrenches_Jacobian(q_7dof, param);
    dJ = dJ_SQUAT_7DOF_Ext_wrenches_Jacobian(q_7dof, dq_7dof, param.L');

    %% calculation of JpQp
    x = zeros(3, param.NbSample);
    dx = zeros(3, param.NbSample);
    ddx = zeros(3, param.NbSample);
    
    switch fctMode
        case 'symbolic'
            x = sym(x);
            dx = sym(dx);
            ddx = sym(ddx);
    end
    
    for ii=1:param.NbSample
        x(:,ii) = Tr(1:3,4,end,ii);
        dx(:,ii) = J([1 2 6],:,ii)*dq_7dof(:,ii); % dx = J dq
        ddx(:,ii) = J([1 2 6],:,ii)*ddq_7dof(:,ii) + dJ([1 2 6],:,ii)*dq_7dof(:,ii); % ddx = J ddq + dJ dq
        
%         COM_local = COM_SQUAT_7DOF(Tr(:, :, :, ii), param);
%         com(:, ii) = COM_local(:, end);
    end
    dddx = calcDeriv(ddx, param.dt_spline);
    
    com = zeros(3, size(q_7dof, 2));
    dcom = zeros(3, size(q_7dof, 2));
    ddcom = zeros(3, size(q_7dof, 2));
    
%     and calculating the forward dynamics 
    [tau, dtau, ddtau, cop, dcop, ddcop] = calc_dynamics_timeseries(q_7dof, dq_7dof, ddq_7dof, param);
    
    
    
    param = update_length(param, 1);
    
    q_1len = sym('q', [param.dof_count 1]);
    dq_1len = sym('dq', [param.dof_count 1]);
%     ddq_1len = sym('ddq', [param.dof_count 1]);
%     dddq_1len = sym('dddq', [param.dof_count 1]);
    
%     lengthQ = 1;
%     q_7dof_1len   = [zeros(param.dofsFromFull(1)-1, lengthQ);   q_1len;   zeros(param.dofsFull(end)-param.dofsFromFull(end), lengthQ)];
%     dq_7dof_1len  = [zeros(param.dofsFromFull(1)-1, lengthQ);   dq_1len;  zeros(param.dofsFull(end)-param.dofsFromFull(end), lengthQ)];
    q_7dof_1len = q_7dof(:, 1);
    dq_7dof_1len = dq_7dof(:, 1);
    
    % generate the inertial and mass matrix
%     H = SQUAT_7DOF_ccg(q_7dof, dq_7dof, param);
    M = SQUAT_7DOF_inm(q_7dof_1len, param); % inertial matrix
    
    % calculate potential energy
    ep = SQUAT_7DOF_ccg(q_7dof, z_7dof, param); % potential energy due to configuration and no velocity
    
    % kinetic energy and geodesic
    switch fctMode
        case 'symbolic'
            dq_3d = reshape(dq_7dof_1len, [1, size(dq_7dof_1len, 1), size(dq_7dof_1len, 2)]);  % vectorize dq properly
            dq_3d = repmat(dq_3d, [size(dq_7dof_1len, 1), 1, 1]);
            
            generic_ek = sum(M.*(dq_3d.^2)); % kinetic energy
            q_dq_total = [q_7dof_1len; dq_7dof_1len];
            
            ek = sym(zeros(size(q_7dof)));
            for ind_x = 1:lengthQ
                q_dq_replace = [q_7dof(:, ind_x); dq_7dof(:, ind_x)];
                ek(:, ind_x) = subs(generic_ek, q_dq_total, q_dq_replace);
            end
            
%             ek = sum(ek_temp_array, 2); % kinetic energy
%             ek = reshape(ek_temp, [size(dq_7dof, 1), size(dq_7dof, 2)]);
            geo = sqrt(ek);
            
        otherwise
            dq_3d = reshape(dq_7dof, [1, size(dq_7dof, 1), size(dq_7dof, 2)]);  % vectorize dq properly
            dq_3d = repmat(dq_3d, [size(dq_7dof, 1), 1, 1]);
            ek_temp = sum(M.*(dq_3d.^2), 2); % kinetic energy
            ek = reshape(ek_temp, [size(dq_7dof, 1), size(dq_7dof, 2)]);
            geo = sqrt(ek);
    end
    
    % 'energy'
    en = dq_7dof .* tau;
    
    % modify the FD so it only keeps the 3 joints we're interested in
    indToUse = param.dofsFromFull;
    tau = tau(indToUse, :);
    dtau = dtau(indToUse, :);
    ddtau = ddtau(indToUse, :);
    ep = ep(indToUse, :);
    ek = ek(indToUse, :);
    geo = geo(indToUse, :);
    en = en(indToUse, :);

    % rearrange all the features into a struct so it's faster to pass
    % around
    feature_out.q = q;
    feature_out.dq = dq;
    feature_out.ddq = ddq;
    feature_out.dddq = dddq;
    feature_out.x = x;
    feature_out.dx = dx;
    feature_out.ddx = ddx;
    feature_out.dddx = dddx;
    feature_out.tau = tau;
    feature_out.dtau = dtau;
    feature_out.ddtau = ddtau;
    feature_out.cop = cop;
    feature_out.dcop = dcop;
    feature_out.ddcop = ddcop;
    feature_out.com = com;
    feature_out.dcom = dcom;
    feature_out.ddcom = ddcom;
    feature_out.ep = ep;
    feature_out.ek = ek;
    feature_out.geo = geo;
    feature_out.en = en;
end