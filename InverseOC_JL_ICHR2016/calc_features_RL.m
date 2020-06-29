function feature_out = calc_features_RL(splineFit, feature_win, param)
    % given the spline variables, produce all the individual variables
    
    calcAngMom = 0;
    
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
    x_max = [];
    
    % Jump-specific variables
    model = inputData(model, q_full(:, 1), dq_full(:, 1), ddq_full(:, 1));
    transformNames = {model.transforms.name};
%     toeTrNum = find(ismember(transformNames,'rAnkle2Toe')==1);
%     toe_pos = [];
%     
%     COM = [];
%     ang_mom = 0;
        
    
      for i = 1:size(q, 2)
        %THIS IS TO SET THEM IN ORDER
        model = inputData(model, q_full(:, i), dq_full(:, i), ddq_full(:, i));

        %MUST CALL THIS TO UPDATE MODEL
        model.forwardPosition();
        model.forwardVelocity();
        model.forwardAcceleration();
        model.inverseDynamics();

        tau_curr = model.torque; %% CF: Tau
        % Cartesian position of end-effector
        t_curr = model.getFrameByName(param.endEffectorName).t;
        x_curr = t_curr(1:3, 4);
        rot_curr = t_curr(1:3, 1:3);
        
        % Distance between end-effector and pick (anchor 1) and place (anchor 2) locations
        x_anchorDiff1(i, :) = norm(x_curr - param.x_anchor_1); %% CF: Distance to target location 1
        x_anchorDiff2(i, :) = norm(x_curr - param.x_anchor_2); %% CF: Distance to target location 2
        
        % Distance between end-effector rotation at pick (anchor 1) and
        % placement (anchor 2) locations
        x_rotDiff1(i, :) = acos((trace(rot_curr*param.rot_anchor_1') -1)/2);
        x_rotDiff2(i, :) = acos((trace(rot_curr*param.rot_anchor_2') -1)/2);
                
        for j = 1:size(param.end_eff_frames, 2)
            ft = model.getFrameByName(param.end_eff_frames{j}).t(1:3, 4);
            % Get Cartesian positions of all anatomical joints (needed for curvature)
            x_listTemp(param.end_eff_framesInds{j}, i) = ft;
            % Compute displacement of end-effector relative to other joints
            x_displace(param.end_eff_framesInds{j}, i) = ft - x_curr; %% CF: Displacement
        end
                
        model.calculateMassMatrix(); % Get inertial matrix for kinecti energy CF         
        M(:, :, i) = model.M;
        
        com(i, :) = com_calc(model); % Get COM position
        
        % CoM and angular momentum
%         if(calcAngMom)
%             if(i==1)
%                 [COM_model_curr, com_body_curr, t_body_curr] = jump2D_com_calc(model);
%                 COM = [COM COM_model_curr];
% 
%             else % approximates velocities between two adjacent positions, assume momentum of first frame is 0
%                 [COM_model_curr, ang_mom_curr, com_body_curr, t_body_curr] = jump2D_com_ang_momentum(model, com_body_prev, t_body_prev);
%                 COM = [COM COM_model_curr];
%                 ang_mom = [ang_mom ang_mom_curr];
%             end
%             com_body_prev = com_body_curr;
%             t_body_prev = t_body_curr;
%             
%         else % onle calc CoM info
%             [COM_model_curr, ~, ~] = jump2D_com_calc(model);
%             COM = [COM COM_model_curr];
% %             ang_mom = [ang_mom 0];
%         end
        
%         model_prev = model; % doesn't work b/c uses pass by reference
        
        
        
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
        
      end
    
%     if(~calcAngMom)
%         ang_mom = zeros(size(COM));
%     end

      % % Getting coordinates of bounding box        
      x_max(1, :) = max(abs(x_listTemp([1 4 7], :)), [], 1);
      x_max(2, :) = max(abs(x_listTemp([2 5 8], :)), [], 1);
      x_max(3, :) = max(abs(x_listTemp([3 6 9], :)), [], 1);

      dx_listTemp = calcDeriv(x_listTemp, param.dt_spline); % Cartesian velocity
      ddx_listTemp = calcDeriv(dx_listTemp, param.dt_spline); % Cartesian acceleration
      dddx_listTemp = calcDeriv(ddx_listTemp, param.dt_spline); % Cartesian jerk
        
        for j = 1:size(param.end_eff_frames, 2)
            x_all(j, :) = normVector(x_listTemp(param.end_eff_framesInds{j}, :)')';    
            dx_all(j, :) = normVector(dx_listTemp(param.end_eff_framesInds{j}, :)')';    
            ddx_all(j, :) = normVector(ddx_listTemp(param.end_eff_framesInds{j}, :)')';
            dddx_all(j, :) = normVector(dddx_listTemp(param.end_eff_framesInds{j}, :)')';
        end
         
        cartCurv = calcCfCartCurvature(dx_listTemp, ddx_listTemp, param.end_eff_framesInds);
        shapeDir = calcCfShapeDirectionSumSqu(x_listTemp, dx_listTemp, ddx_listTemp, param.end_eff_framesInds);
        
        dx_displace = calcDeriv(x_displace, param.dt_spline);
        ddx_displace = calcDeriv(dx_displace, param.dt_spline);
        dddx_displace = calcDeriv(ddx_displace, param.dt_spline);
        
        for j = 1:size(param.end_eff_frames, 2)
            x_displaceNorm(j, :) = normVector(x_displace(param.end_eff_framesInds{j}, :)')';
            dx_displaceNorm(j, :) = normVector(dx_displace(param.end_eff_framesInds{j}, :)')';
            ddx_displaceNorm(j, :) = normVector(ddx_displace(param.end_eff_framesInds{j}, :)')';
            dddx_displaceNorm(j, :) = normVector(dddx_displace(param.end_eff_framesInds{j}, :)')';
        end
    
    dx = calcDeriv(x, param.dt_spline); % CF: End-effector velocity
    ddx = calcDeriv(dx, param.dt_spline); % CF: End-effector acceleration
    dddx = calcDeriv(ddx, param.dt_spline); % CF: End-effector jerk
    x_cartCurv = calcCfCartCurvature(dx, ddx, {1:3});
        
    dtau = calcDeriv(tau, param.dt_spline); % CF: Torque change
    ddtau = calcDeriv(dtau, param.dt_spline); % CF: Torque effort
%     dcop = calcDeriv(cop, param.dt_spline);
%     ddcop = calcDeriv(dcop, param.dt_spline);
    dcom = calcDeriv(com, param.dt_spline);
    ddcom = calcDeriv(dcom, param.dt_spline);
    
%     dCOM = calcDeriv(COM, param.dt_spline); % COM cartesian trajectories
%     ddCOM = calcDeriv(dCOM, param.dt_spline);
%     dddCOM = calcDeriv(ddCOM, param.dt_spline);
%     dCOM_mag = sqrt(sum(dCOM.^2, 1)); % COM magnitude trajectories
%     ddCOM_mag = sqrt(sum(ddCOM.^2, 1));
%     dddCOM_mag = sqrt(sum(dddCOM.^2, 1));
    
%     ang_mom(1) = ang_mom(2); % reset first value to 2nd value so derivatives aren't massive over first frame
%     dang_mom = calcDeriv(ang_mom, param.dt_spline);
%     ddang_mom = calcDeriv(dang_mom, param.dt_spline);
    
    dq_3dof = reshape(dq_full, [1, size(dq_full, 1), size(dq_full, 2)]);  % vectorize dq properly
    dq_3dof = repmat(dq_3dof, [size(dq_full, 1), 1, 1]);
    
    % CF: Kinetic Energy
    ek_temp = sum(M.*(dq_3dof.^2), 2); % kinetic energy
    ek = reshape(ek_temp, [size(dq_full, 1), size(dq_full, 2)]);
    % CF: Squared root Kinetic Energy
    geo = sqrt(ek);
    
    % CF: Power
    en = dq_full .* tau;
          
    indToUse = param.dofsFromFull;
    tau = tau(indToUse, :);
    dtau = dtau(indToUse, :);
    ddtau = ddtau(indToUse, :);
%     ep = ep(indToUse, :);
    ek = ek(indToUse, :);
    geo = geo(indToUse, :);
    en = en(indToUse, :);  
    
%     dofHardcode = indToUse;
%     ddq_tp1 = zeros(length(dofHardcode), lengthQ);
%     dq_tp1 = zeros(length(dofHardcode), lengthQ);
%     q_tp1 = zeros(length(dofHardcode), lengthQ);
    
    % All global X pos features should be relative to target/toe
%     dToeX = calcDeriv(toe_pos(1,:), param.dt_spline); 
%     ddToeX = calcDeriv(dToeX, param.dt_spline);
%     dToeZ = calcDeriv(toe_pos(3,:), param.dt_spline);
%     ddToeZ = calcDeriv(dToeZ, param.dt_spline);
    
    % Distances in forward jump direction (Global X)
%     ToeTargX = abs(toe_pos(1,:) - param.jump.locationLand);
%     dToeTargX = calcDeriv(ToeTargX, param.dt_spline);
%     COMToeX = abs(COM(1,:) - toe_pos(1,:));
%     dCOMToeX = calcDeriv(COMToeX, param.dt_spline);
%     COMTargX = abs(COM(1,:) - param.jump.locationLand);
%     dCOMTargX = calcDeriv(COMTargX, param.dt_spline);
    
    
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
    
    feature_out.x_anchor_1 = x_anchorDiff1;
    feature_out.x_anchor_2 = x_anchorDiff2;
    feature_out.rot_anchor_1 = x_rotDiff1;
    feature_out.rot_anchor_2 = x_rotDiff2;
    
    feature_out.x_max = x_max;
    feature_out.x_all = x_all;
    feature_out.dx_all = dx_all;
    feature_out.ddx_all = ddx_all;
    feature_out.dddx_all = dddx_all;
    
    feature_out.x_displace = x_displaceNorm;
    feature_out.dx_displace = dx_displaceNorm;
    feature_out.ddx_displace = ddx_displaceNorm;
    feature_out.dddx_displace = dddx_displaceNorm;
    
    feature_out.cartCurv = cartCurv;
    feature_out.shapeDir = shapeDir;
    feature_out.x_cartCurv = x_cartCurv;
    
    
    feature_out.tau = tau; % (close)
    feature_out.dtau = dtau;
    feature_out.ddtau = ddtau;
%     feature_out.cop = cop;
%     feature_out.dcop = dcop;
%     feature_out.ddcop = ddcop;
    feature_out.com = com;
    feature_out.dcom = dcom;
    feature_out.ddcom = ddcom;
%     feature_out.ep = ep;
    feature_out.ek = ek;
    feature_out.geo = geo; % coresponds to "kinetic_en_sqrt"
    feature_out.en = en;
    
%     feature_out.ddq_tp1 = ddq_tp1;
%     feature_out.dq_tp1 = dq_tp1;
%     feature_out.q_tp1 = q_tp1;
    
%     feature_out.COM = COM;
%     feature_out.dCOM = dCOM;
%     feature_out.ddCOM = ddCOM;
%     feature_out.dddCOM = dddCOM;
%     feature_out.dCOM_mag = dCOM_mag;
%     feature_out.ddCOM_mag = ddCOM_mag;
%     feature_out.dddCOM_mag = dddCOM_mag;
%     feature_out.COM_height = COM(3,:);
%     feature_out.dCOM_height = dCOM(3,:);
%     feature_out.ddCOM_height = ddCOM(3,:);
%     feature_out.dddCOM_height = dddCOM(3,:);
    
%     feature_out.dToeZ = dToeZ;
%     feature_out.ddToeZ = ddToeZ;
    
%     feature_out.ToeTargX = ToeTargX;
%     feature_out.dToeTargX = dToeTargX;
%     feature_out.COMToeX = COMToeX;
%     feature_out.dCOMToeX = dCOMToeX;
%     feature_out.COMTargX = COMTargX;
%     feature_out.dCOMTargX = dCOMTargX;
    
    % ddq, dddq, tau, dq_tau (en), kinetic_en_sqrt (geo) 
    % divided for torso, arms, legs
%     feature_out.ddq_tor = ddq(1:4,:);
%     feature_out.dddq_tor = dddq(1:4,:);
%     feature_out.tau_tor = tau(1:4,:);
%     feature_out.en_tor = en(1:4,:);
%     feature_out.geo_tor = geo(1:4,:);
%     
%     feature_out.ddq_arm = ddq(5:8,:);
%     feature_out.dddq_arm = dddq(5:8,:);
%     feature_out.tau_arm = tau(5:8,:);
%     feature_out.en_arm = en(5:8,:);
%     feature_out.geo_arm = geo(5:8,:);
    
%     feature_out.ddq_leg = ddq(9:14,:);
%     feature_out.dddq_leg = dddq(9:14,:);
%     feature_out.tau_leg = tau(9:14,:);
%     feature_out.en_leg = en(9:14,:);
%     feature_out.geo_leg = geo(9:14,:);
    
    % centroidal angular momentum
%     feature_out.ang_mom = ang_mom;
%     feature_out.dang_mom = dang_mom;
%     feature_out.ddang_mom = ddang_mom;  
end

function cf = calcCfCartCurvature(dx, ddx, frameInds)
    lenTime = size(dx, 2);
    lenDofs = length(frameInds);

    cf = zeros(lenDofs, lenTime);
    for i = 1:lenTime
        for j = 1:lenDofs
            inds = frameInds{j};
            num = norm(cross(ddx(inds, i), dx(inds, i))); % todo can vectorize
            den = norm(dx(inds, i).^(3/2));
            if den > 0
                cfVal = num/den; % if there is no dx, then we're dividing by zero
            else
                cfVal = 0;
            end

            cf(j, i) = cfVal;
        end
    end
end

function cf = calcCfShapeDirectionSumSqu(x, dx, ddx, frameInds)

cf = 0;
% return;

lenTime = size(dx, 2);
lenDofs = length(frameInds);

            cf = zeros(lenDofs, lenTime);
            for i = 2:lenTime
                for j = 1:lenDofs
                    % assemble
                    inds = frameInds{j};
                    tempx(j, :) = x(inds, i);
                    tempdx(j, :) = dx(inds, i);
                    tempddx(j, :) = ddx(inds, i);
                end
                
                coeff = pca(tempx)';
                
                pcadx = zeros(lenDofs, 2);
                pcaddx = zeros(lenDofs, 2);
                for j = 1:lenDofs
                    pcadx(j, :) = coeff*tempdx(j, :)';
                    pcaddx(j, :) = coeff*tempddx(j, :)';
                end
           
                num = (pcaddx(2)*pcadx(1) - pcaddx(1)*pcadx(2)) .^ 2;
                den = ((pcadx(1) + pcadx(2)) .^2) .^ (3/2);
                
                if den ~= 0
                    cf(i, j) = num/den;
                else
                    cf(i, j) = 0;
                end
            end
end
