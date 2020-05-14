function feature_out = calc_features(splineFit, feature_win, param)
    % given the spline variables, produce all the individual variables
        
    % calculate the joint angles
    if isstruct(splineFit) 
        if isfield(splineFit.q, 'sp_qd0') % spline fit is actually a spline
            [q, dq, ddq, dddq] = calc_q(splineFit.q, param);
            
            fctMode = 'spline';
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
    
    efToUse = param.dofsFromFull(end)+2; % +1 to get the next EF entry, then another +1 since there's a gap EF in the definition
    
    for ii=1:param.NbSample
        x(:,ii) = Tr(1:3,4,efToUse,ii);
        dx(:,ii) = J([1 2 6],:,ii)*dq_7dof(:,ii); % dx = J dq
        ddx(:,ii) = J([1 2 6],:,ii)*ddq_7dof(:,ii) + dJ([1 2 6],:,ii)*dq_7dof(:,ii); % ddx = J ddq + dJ dq
        
%         COM_local = COM_SQUAT_7DOF(Tr(:, :, :, ii), param);
%         com(:, ii) = COM_local(:, end);
    end
    dddx = calcDeriv(ddx, param.dt_spline);
    
    com = zeros(3, size(q_7dof, 2));
    dcom = calcDeriv(com, param.dt_spline);
    ddcom = calcDeriv(dcom, param.dt_spline);
    
%     and calculating the forward dynamics 
    [tau, dtau, ddtau, cop, dcop, ddcop] = calc_dynamics_timeseries(q_7dof, dq_7dof, ddq_7dof, param);
    
    % generate the inertial and mass matrix
%     H = SQUAT_7DOF_ccg(q_7dof, dq_7dof, param);
    M = SQUAT_7DOF_inm(q_7dof, param); % inertial matrix
    
    % calculate potential energy
    ep = SQUAT_7DOF_ccg(q_7dof, z_7dof, param.base0, param); % potential energy due to configuration and no velocity
    
    % kinetic energy and geodesic
    switch fctMode
        case 'symbolic'
            ek = zeros(size(q_7dof));
            geo = zeros(size(q_7dof));
            
        otherwise
            dq_3d = reshape(dq_7dof, [1, size(dq_7dof, 1), size(dq_7dof, 2)]);  % vectorize dq properly
            dq_3d = repmat(dq_3d, [size(dq_7dof, 1), 1, 1]);
            ek_temp = sum(M.*(dq_3d.^2), 2); % kinetic energy
            ek = reshape(ek_temp, [size(dq_7dof, 1), size(dq_7dof, 2)]);
            geo = sqrt(ek);
    end
    
    % 'energy'
    en = dq_7dof .* tau;
    
    indToUse = param.dofsFromFull;
    
    % constraints for the next timestep
    tau_delta = tau - SQUAT_7DOF_dyn3(q_7dof, dq_7dof, z_7dof, param.base0, param.ext_wrenches, param);
%     dofHardcode = 1:4;
    dofHardcode = indToUse;
    ddq_tp1 = zeros(length(dofHardcode), lengthQ);
    dq_tp1 = zeros(length(dofHardcode), lengthQ);
    q_tp1 = zeros(length(dofHardcode), lengthQ);
    for ind_const = 1:lengthQ
        invM = inv(M(dofHardcode, dofHardcode, ind_const));
%         invM = eye(4);
        
        ddq_tp1(:, ind_const) = invM*(tau_delta(dofHardcode, ind_const));
        dq_tp1(:, ind_const)  = dq_7dof(dofHardcode, ind_const) + param.dt_full*ddq_7dof(dofHardcode, ind_const);
        q_tp1(:, ind_const)   =  q_7dof(dofHardcode, ind_const) + param.dt_full* dq_7dof(dofHardcode, ind_const);
    end
    
    % modify the FD so it only keeps the 3 joints we're interested in

    tau = tau(indToUse, :);
    dtau = dtau(indToUse, :);
    ddtau = ddtau(indToUse, :);
    ep = ep(indToUse, :);
    ek = ek(indToUse, :);
    geo = geo(indToUse, :);
    en = en(indToUse, :);   
    
%     ddq_tp1 = ddq_tp1(indToUse, :);
%     dq_tp1 = dq_tp1(indToUse, :);
%     q_tp1 = q_tp1(indToUse, :);
    
%     % TESTING
%     tau_test = SQUAT_7DOF_dyn3(q_7dof, dq_7dof, ddq_7dof, param.base0, param.ext_wrenches, param);
%     M_test = SQUAT_7DOF_inm(q_7dof, param);
%     ddq_test = ddq_7dof;
%     ccg_test = SQUAT_7DOF_ccg(q_7dof, dq_7dof, param.base0, param);
%     
% %     for indlala = 1:201
% %         mddq_array(:, indlala) = M_test(:, :, indlala) * ddq_test(:, indlala);
% %     end
%     
%     ddq_3d = reshape(ddq_test, [1, size(ddq_test, 1), size(ddq_test, 2)]);  % vectorize dq properly
%     ddq_3d = repmat(ddq_3d, [size(ddq_test, 1), 1, 1]);
%     Mddq_test = sum(M_test.*(ddq_3d), 2); % kinetic energy
%     Mddq2_test = reshape(Mddq_test, [size(ddq_test, 1), size(ddq_test, 2)]);
%     tau2_test = Mddq2_test + ccg_test;
    
%     figure; plot(tau_test'); hold on; plot(tau2_test', '.')
%     figure; plot(ddq'); hold on; plot(ddq_tp1', '.');
    
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
    
    feature_out.ddq_tp1 = ddq_tp1;
    feature_out.dq_tp1 = dq_tp1;
    feature_out.q_tp1 = q_tp1;    
end