function feature_out = unnormalize_features(feature_use, param)

    keyboard;

    % normalize
    feature_out.q   =  feature_use.q   ./repmat(param.coeff_q,     1, size(feature_use.q, 2));
    feature_out.dq  =  feature_use.dq  ./repmat(param.coeff_dq,    1, size(feature_use.dq, 2));
    feature_out.ddq =  feature_use.ddq ./repmat(param.coeff_ddq,   1, size(feature_use.ddq, 2));
    feature_out.dddq = feature_use.dddq./repmat(param.coeff_dddq,  1, size(feature_use.dddq, 2));
    
    feature_out.x   =  feature_use.x   ./repmat(param.coeff_x,     1, size(feature_use.x, 2));
    feature_out.dx  =  feature_use.dx  ./repmat(param.coeff_dx,    1, size(feature_use.dx, 2));
    feature_out.ddx =  feature_use.ddx ./repmat(param.coeff_ddx,   1, size(feature_use.ddx, 2));
    feature_out.dddx = feature_use.dddx./repmat(param.coeff_dddx,  1, size(feature_use.dddx, 2));
    
    feature_out.tau =  feature_use.tau ./repmat(param.coeff_tau,    1, size(feature_use.tau, 2));
    feature_out.dtau = feature_use.dtau./repmat(param.coeff_dtau,   1, size(feature_use.dtau, 2));
end