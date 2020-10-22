function feature_opt = feature_windowed(feature_full, ind1, ind2, param)
    % given the collated dataset, we want to crop and resample data as
    % needed    
%     param.windowMode = 'wide';
%     
    switch  param.windowMode 
        case 'narrow'
            t_windowCount =   feature_full.t(:, ind1:ind2);
            q_windowCount =   feature_full.q(:, ind1:ind2);
            dq_windowCount =  feature_full.dq(:, ind1:ind2);
            ddq_windowCount = feature_full.ddq(:, ind1:ind2);
            dddq_windowCount = feature_full.dddq(:, ind1:ind2);
            
            % resample to ensure the same length as the spline length
%             dt_opt = (t_windowCount(end)-t_windowCount(1))/(param.spline_length-1);   
            dt_opt = param.dt_spline;
            t_opt = t_windowCount(1):dt_opt:(t_windowCount(end));
            t   = t_opt;
            
            if sqrt(sum((t_windowCount - t_opt).^2)) > 1e-3 % check to see if the array is actually different via RMS
                feature_narrow.q   = interp1(t_windowCount', q_windowCount', t_opt')';
                feature_narrow.dq   = interp1(t_windowCount', dq_windowCount', t_opt')';
                feature_narrow.ddq   = interp1(t_windowCount', ddq_windowCount', t_opt')';
                feature_narrow.dddq   = interp1(t_windowCount', dddq_windowCount', t_opt')';
            else
                feature_narrow.q   = q_windowCount;
                feature_narrow.dq   = dq_windowCount;
                feature_narrow.ddq   = ddq_windowCount;
                feature_narrow.dddq   = dddq_windowCount;
            end

            % now calculate the remaining features
            feature_opt = calc_features(feature_narrow, [], param);
            feature_opt.t = t;

            %     feature_opt.dq  = interp1(t_windowCount', dq_windowCount', t_opt')';
            %     feature_opt.ddq = interp1(t_windowCount', ddq_windowCount', t_opt')';
            
        case 'wide'
            offsetVal = 3;
            ind_narrow = ind1:ind2;
            ind_wide = (ind1-offsetVal):(ind2+offsetVal);
            ind_wide2narrow = (offsetVal+1):(length(ind_wide)-offsetVal);
%             ind_narrow
%             ind_wide(ind_wide2narrow)
%             ind_wide = [ind1 ind1 ind1 ind1:ind2 ind2 ind2 ind2];
            len_t = length(feature_full.t);
            
            ind_wide(ind_wide < 1) = 1;
            ind_wide(ind_wide > len_t) = len_t;
            
%             param_wide = update_length(length(ind_wide), length(ind_wide), param);

            % pass it in for the wider window for now
            feature_wide.q =   feature_full.q(:, ind_wide);
            feature_wide.dq =  feature_full.dq(:, ind_wide);
            feature_wide.ddq = feature_full.ddq(:, ind_wide);
            feature_wide.dddq = feature_full.dddq(:, ind_wide);
            feature_wide = calc_features(feature_wide, [], param); % it will come back out the proper length
            
            % crop down to narrow length
            feature_narrow = feature_crop(feature_wide, ind_wide2narrow);

            % resample to ensure the same length as the spline length
            t_windowCount =   feature_full.t(:, ind_narrow);
            
            dt_opt = (t_windowCount(end)-t_windowCount(1))/(param.spline_length-1);
            t_opt = t_windowCount(1):dt_opt:(t_windowCount(end));
            
            % perform interpolation and crop it
            if sqrt(sum((t_windowCount - t_opt).^2)) > 1e-3 % check to see if the array is actually different via RMS
                % no interpolation needed
                feature_opt = feature_interpolation(feature_narrow, t_windowCount, t_opt);
            else
                feature_opt = feature_narrow;
            end
            feature_opt.t = t_opt;
    end
end

function feature_out = feature_crop(feature_wide, ind_out)
    feature_out.q  =   feature_wide.q(:, ind_out);
    feature_out.dq =   feature_wide.dq(:, ind_out);
    feature_out.ddq =  feature_wide.ddq(:, ind_out);
    feature_out.dddq = feature_wide.dddq(:, ind_out);
%     feature_out.x =    feature_wide.x(:, ind_out);
%     feature_out.dx =   feature_wide.dx(:, ind_out);
%     feature_out.ddx =  feature_wide.ddx(:, ind_out);
%     feature_out.dddx = feature_wide.dddx(:, ind_out);
%     feature_out.tau =  feature_wide.tau(:, ind_out);
%     feature_out.dtau = feature_wide.dtau(:, ind_out);
%     feature_out.ddtau =feature_wide.ddtau(:, ind_out);
%     feature_out.cop =  feature_wide.cop(:, ind_out);
%     feature_out.dcop = feature_wide.dcop(:, ind_out);
%     feature_out.ddcop =feature_wide.ddcop(:, ind_out);
%     feature_out.com =  feature_wide.com(:, ind_out);
%     feature_out.dcom = feature_wide.dcom(:, ind_out);
%     feature_out.ddcom =feature_wide.ddcom(:, ind_out);
%     feature_out.ep =   feature_wide.ep(:, ind_out);
%     feature_out.ek =   feature_wide.ek(:, ind_out);
%     feature_out.geo =  feature_wide.geo(:, ind_out);
%     feature_out.en =   feature_wide.en(:, ind_out);

    feature_out.joint_limit = feature_wide.joint_limit(:, ind_out);
%     feature_out.manip_rot = feature_wide.manip_rot(:, ind_out);
%     feature_out.manip_trans = feature_wide.manip_trans(:, ind_out);
    feature_out.manipulability = feature_wide.manipulability(:, ind_out);
%     feature_out.half_joint_task = feature_wide.half_joint_task(:, ind_out);
    feature_out.orientation_length = feature_wide.orientation_length(:, ind_out);
    feature_out.joint_length = feature_wide.joint_length(:, ind_out);
    feature_out.task_length = feature_wide.task_length(:, ind_out);
    feature_out.task_jerk = feature_wide.task_jerk(:, ind_out);
    feature_out.joint_jerk = feature_wide.joint_jerk(:, ind_out);
end

function feature_out = feature_interpolation(feature_wide, t_orig, t_opt)
    feature_out.q  =   interp1(t_orig', feature_wide.q', t_opt')';
    feature_out.dq =   interp1(t_orig', feature_wide.dq', t_opt')';
    feature_out.ddq =  interp1(t_orig', feature_wide.ddq', t_opt')';
    feature_out.dddq = interp1(t_orig', feature_wide.dddq', t_opt')';
%     feature_out.x =    interp1(t_orig', feature_wide.x', t_opt')';
%     feature_out.dx =   interp1(t_orig', feature_wide.dx', t_opt')';
%     feature_out.ddx =  interp1(t_orig', feature_wide.ddx', t_opt')';
%     feature_out.dddx = interp1(t_orig', feature_wide.dddx', t_opt')';
%     feature_out.tau =  interp1(t_orig', feature_wide.tau', t_opt')';
%     feature_out.dtau = interp1(t_orig', feature_wide.dtau', t_opt')';
%     feature_out.ddtau =interp1(t_orig', feature_wide.ddtau', t_opt')';
%     feature_out.cop =  interp1(t_orig', feature_wide.cop', t_opt')';
%     feature_out.dcop = interp1(t_orig', feature_wide.dcop', t_opt')';
%     feature_out.ddcop =interp1(t_orig', feature_wide.ddcop', t_opt')';
%     feature_out.com =  interp1(t_orig', feature_wide.com', t_opt')';
%     feature_out.dcom = interp1(t_orig', feature_wide.dcom', t_opt')';
%     feature_out.ddcom =interp1(t_orig', feature_wide.ddcom', t_opt')';
%     feature_out.ep =   interp1(t_orig', feature_wide.ep', t_opt')';
%     feature_out.ek =   interp1(t_orig', feature_wide.ek', t_opt')';
%     feature_out.geo =  interp1(t_orig', feature_wide.geo', t_opt')';
%     feature_out.en =   interp1(t_orig', feature_wide.en', t_opt')';
end