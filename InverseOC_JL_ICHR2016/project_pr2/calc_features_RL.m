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
    
    mdl = param.model;
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
    
    for i = 1:size(q, 2)
        % updating cost functions
        joint_limit(i, :) = joint_limit_objective(mdl, q_full, dq_full, ddq_full, i);
%         manip_rot(i, :) = manip_rot_objective(mdl, q_full, dq_full, ddq_full, i);
%         manip_trans(i, :) = manip_trans_objective(mdl, q_full, dq_full, ddq_full, i);
        manipulability(i, :) = manipulability_objective(mdl, q_full, dq_full, ddq_full, i);
%         task_jerk(i,:) = task_jerk_objective(mdl, q_full, dq_full, ddq_full, i);
%         joint_jerk(i,:) = joint_jerk_objective(mdl, q_full, dq_full, ddq_full, i);
        
        if i > 1
%             half_joint_task(i, :) = half_joint_task_objective(mdl, q_full, dq_full, ddq_full, i-1, i);
            orientation_length(i, :) = orientation_length_objective(mdl, q_full, dq_full, ddq_full, i-1, i);
            joint_length(i, :) = joint_length_objective(mdl, q_full, dq_full, ddq_full, i-1, i);
            task_length(i, :) = task_length_objective(mdl, q_full, dq_full, ddq_full, i-1, i);
            
        end
    end
    task_jerk = task_jerk_objective(mdl, q_full, dq_full, ddq_full);
    % rearrange all the features into a struct so it's faster to pass
    % around
    feature_out.q = q;
    feature_out.dq = dq;
    feature_out.ddq = ddq;
    feature_out.dddq = dddq;
   
    feature_out.joint_limit = joint_limit';
%     feature_out.manip_rot = manip_rot';
%     feature_out.manip_trans = manip_trans';
    feature_out.manipulability = manipulability';
%     feature_out.half_joint_task = half_joint_task';
    feature_out.orientation_length = orientation_length';
    feature_out.joint_length = joint_length';
    feature_out.task_length = task_length';
    feature_out.task_jerk = task_jerk;
    feature_out.joint_jerk = dddq;
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
