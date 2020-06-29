function curr_spline_fit = applyKnotOffsetAndSpline(curr_q_knot, feature_opt, ii_vec, ii_h, param)
    curr_q_knot(ii_vec) = curr_q_knot(ii_vec) + param.h_array(ii_h);
    
    y_opt_mat = vec2mat_multi(curr_q_knot, param.matVec_struct);
    curr_spline_fit =  calc_direct_spline(y_opt_mat, param);
    
% % %     % dq/ddq end condition constraints should proprogate through from the
% % %     % knot shift (h_array) here, so this code runs the spline calculation
% % %     % twice to try to account for it. 
% % %     curr_q_knot(ii_dofs, ii_knot) = curr_q_knot(ii_dofs, ii_knot) + param.h_array(ii_h);
% % %     old_spline_fit.q =  calc_direct_spline(curr_q_knot, param);
% % %     
% % %     feature_use = calc_features(old_spline_fit, [], param);
% % %     feature_use.t = feature_opt.t;
% % %     param = set_constraints(feature_use, param);
% % %     curr_spline_fit =  calc_direct_spline(curr_q_knot, param);
% % %     
%     max(max(abs(old_spline_fit.q.sp_qd0.coefs -  curr_spline_fit.sp_qd0.coefs)))
end