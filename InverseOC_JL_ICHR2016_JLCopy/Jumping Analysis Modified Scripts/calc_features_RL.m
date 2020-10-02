function feature_out = calc_features_RL(splineFit, feature_win, param)
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
    q_full = zeros(length(param.dofsFull), lengthQ);
    q_full(param.dofsFromFull, :) = q;
    dq_full = zeros(length(param.dofsFull), lengthQ);
    dq_full(param.dofsFromFull, :) = dq;
    ddq_full = zeros(length(param.dofsFull), lengthQ);
    ddq_full(param.dofsFromFull, :) = ddq;
    
%     q_full   = [zeros(param.dofsFromFull(1)-1, lengthQ);   q;   zeros(param.dofsFull(end)-param.dofsFromFull(end), lengthQ)];
%     dq_full  = [zeros(param.dofsFromFull(1)-1, lengthQ);   dq;  zeros(param.dofsFull(end)-param.dofsFromFull(end), lengthQ)];
%     ddq_full = [zeros(param.dofsFromFull(1)-1, lengthQ);   ddq; zeros(param.dofsFull(end)-param.dofsFromFull(end), lengthQ)];
    z_3dof = zeros(3, size(q, 2));     
    z_full = zeros(size(q_full));
    z_vec = zeros(size(q_full, 1), 1);
    
    model = param.model;
    x = [];
    dx = [];
    ddx = [];
    dddx = [];
    tau = [];
    dtau = [];
    ddtau = [];
    cop = [];
    dcop = [];
    ddcop = [];
    com = [];
    dcom = [];
    ddcom = [];
    ep = [];
    ek = [];
    geo = [];
    en = [];
    
    % Jump-specific variables
    model = inputData(model, q_full(:, 1), dq_full(:, 1), ddq_full(:, 1));
    transformNames = {model.transforms.name};
    toeTrNum = find(ismember(transformNames,'rAnkle2Toe')==1);
    toe_pos = [];
    
    COM = [];
    
    
      for i = 1:size(q, 2)
        %THIS IS TO SET THEM IN ORDER
        model = inputData(model, q_full(:, i), dq_full(:, i), ddq_full(:, i));

        %MUST CALL THIS TO UPDATE MODEL
        model.forwardPosition();
        model.forwardVelocity();
        model.forwardAcceleration();
        model.inverseDynamics();

        tau_curr = model.torque;
        x_curr = model.getFrameByName(param.endEffectorName).t(1:3, 4);
%   x_curr = [0 0 0]';
        
        % Jump-specific variables
        toe_pos_curr = model.transforms(toeTrNum).frame_out.t(1:3,4);
        COM_curr = jump2D_com_calc(model);

        model.calculateMassMatrix();          %  model = inputData(model, q(:, i)', dq(:, i)', ddq(:, i)');
        M(:, :, i) = model.M;
        
%         model.calculateCentrifugalCoriolis();   model = inputData(model, q(:, i)', z_vec, z_vec);
%         Cdq = model.V;
%         
%         model.calculateGravity();               model = inputData(model, q(:, i)', z_vec, z_vec);
%         G = model.G;
%     
%         ep_curr = Cdq + G;
    
        x =     [x      x_curr];
%         dx =    [dx     dx_curr];
%         ddx =   [ddx    ddx_curr];
%         dddx =  [ddddx  dddx_curr];
        tau =   [tau    tau_curr];
%         dtau =  [dtau   dtau_curr];
%         ddtau = [ddtau  ddtau_curr];
%         cop =   [cop    zeros3];
%         dcop =  [dcop   ];
%         ddcop = [ddcop  ];
%         com =   [com    zeros3];
%         dcom =  [dcom ];
%         ddcom = [ddcom];
%         ep =    [ep     ep_curr];
%         ek =    [ek];
%         geo =   [geo];
%         en =    [en];
        toe_pos = [toe_pos toe_pos_curr];
        
        COM = [COM COM_curr];
        
      end
    
      cop = z_3dof;
      com = z_3dof;
    
    dx = calcDeriv(x, param.dt_spline);
    ddx = calcDeriv(dx, param.dt_spline);
    dddx = calcDeriv(ddx, param.dt_spline);
    dtau = calcDeriv(tau, param.dt_spline);
    ddtau = calcDeriv(dtau, param.dt_spline);
    dcop = calcDeriv(cop, param.dt_spline);
    ddcop = calcDeriv(dcop, param.dt_spline);
    dcom = calcDeriv(com, param.dt_spline);
    ddcom = calcDeriv(dcom, param.dt_spline);
    
    dCOM = calcDeriv(COM, param.dt_spline); % COM cartesian trajectories
    ddCOM = calcDeriv(dCOM, param.dt_spline);
    dddCOM = calcDeriv(ddCOM, param.dt_spline);
    dCOM_mag = sqrt(sum(dCOM.^2, 1)); % COM magnitude trajectories
    ddCOM_mag = sqrt(sum(ddCOM.^2, 1));
    dddCOM_mag = sqrt(sum(dddCOM.^2, 1));
    
    
    dq_3dof = reshape(dq_full, [1, size(dq_full, 1), size(dq_full, 2)]);  % vectorize dq properly
    dq_3dof = repmat(dq_3dof, [size(dq_full, 1), 1, 1]);
    ek_temp = sum(M.*(dq_3dof.^2), 2); % kinetic energy
    ek = reshape(ek_temp, [size(dq_full, 1), size(dq_full, 2)]);
    geo = sqrt(ek);
    en = dq_full .* tau;
      
    ep = z_full;
    
    indToUse = param.dofsFromFull;
    tau = tau(indToUse, :);
    dtau = dtau(indToUse, :);
    ddtau = ddtau(indToUse, :);
    ep = ep(indToUse, :);
    ek = ek(indToUse, :);
    geo = geo(indToUse, :);
    en = en(indToUse, :);  
    
    dofHardcode = indToUse;
    ddq_tp1 = zeros(length(dofHardcode), lengthQ);
    dq_tp1 = zeros(length(dofHardcode), lengthQ);
    q_tp1 = zeros(length(dofHardcode), lengthQ);
    
    % All global X pos features should be relative to target/toe
%     dToeX = calcDeriv(toe_pos(1,:), param.dt_spline); 
%     ddToeX = calcDeriv(dToeX, param.dt_spline);
    dToeZ = calcDeriv(toe_pos(3,:), param.dt_spline);
    ddToeZ = calcDeriv(dToeZ, param.dt_spline);
    
    % Distances in forward jump direction (Global X)
    ToeTargX = abs(toe_pos(1,:) - param.jump.locationLand);
    dToeTargX = calcDeriv(ToeTargX, param.dt_spline);
    COMToeX = abs(COM(1,:) - toe_pos(1,:));
    dCOMToeX = calcDeriv(COMToeX, param.dt_spline);
    COMTargX = abs(COM(1,:) - param.jump.locationLand);
    dCOMTargX = calcDeriv(COMTargX, param.dt_spline);
    
    
%     % CoM velocity at landing (projectile motion)
%     com_land_slope = (com_land_post(3) - com_land(3)) / (com_land_post(1) - com_land(1));
%     com_land_offset = com_land(3) - com_land_slope*com_land(1);
%     com_land_predictedX = com_land_offset/com_land_slope;
%     cost_com_land_angle_to_targ = sin(com_land(3)/(param.jump.locationLand - com_land(1)));
    
    
    % rearrange all the features into a struct so it's faster to pass
    % around
    feature_out.q = q;
    feature_out.dq = dq;
    feature_out.ddq = ddq;
    feature_out.dddq = dddq;
    feature_out.x = x;
    feature_out.dx = dx;
    feature_out.ddx = ddx; % ISSUE
    feature_out.dddx = dddx; % ISSUE
    feature_out.tau = tau; % (close)
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
    
    feature_out.COM = COM;
    feature_out.dCOM = dCOM;
    feature_out.ddCOM = ddCOM;
    feature_out.dddCOM = dddCOM;
    feature_out.dCOM_mag = dCOM_mag;
    feature_out.ddCOM_mag = ddCOM_mag;
    feature_out.dddCOM_mag = dddCOM_mag;
    feature_out.COM_height = COM(3,:);
    feature_out.dCOM_height = dCOM(3,:);
    feature_out.ddCOM_height = ddCOM(3,:);
    feature_out.dddCOM_height = dddCOM(3,:);
    
    feature_out.dToeZ = dToeZ;
    feature_out.ddToeZ = ddToeZ;
    
    feature_out.ToeTargX = ToeTargX;
    feature_out.dToeTargX = dToeTargX;
    feature_out.COMToeX = COMToeX;
    feature_out.dCOMToeX = dCOMToeX;
    feature_out.COMTargX = COMTargX;
    feature_out.dCOMTargX = dCOMTargX;
    
end

