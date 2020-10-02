function J = calc_diff_num(x, param, diffType, calcMode, extraParam)
%     diffType = 'symmetrical'; % forward backward central symmetrical
    
    switch diffType
        case 'forward'
            h_offset1 = +1;
            h_offset2 =  0;
            
        case 'backward'
            h_offset1 =  0;
            h_offset2 = -1;
            
        case 'central'
            h_offset1 = +0.5;
            h_offset2 = -0.5;
            
        case 'symmetrical'
            h_offset1 = +1;
            h_offset2 = -1;
    end
    
    h_ind1 = param.h_vals == h_offset1;
    h_ind2 = param.h_vals == h_offset2;
    
    switch calcMode
        case 'cost'
            x_ph_fct = x(:, :, h_ind1); % expects the cost function itself
            x_mh_fct = x(:, :, h_ind2);
            
        case 'const'
            x_ph_fct = x(extraParam.dofVal, extraParam.knotVal, h_ind1) - extraParam.y_base;
            x_mh_fct = x(extraParam.dofVal, extraParam.knotVal, h_ind2) - extraParam.y_base;
    end
    
    len = size(x_ph_fct, 2);
    num = (x_ph_fct - x_mh_fct);
    den = param.h_array(h_ind1) - param.h_array(h_ind2);
    num_diff = num / den; % sum of (ddq(x+h) - ddq(x-h))/2h
    J = sum(num_diff, 2) / len; % summing over all units of time, divided by time
end