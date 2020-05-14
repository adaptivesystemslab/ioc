function rmse_vals = calc_rmse(rmse_fct, feature_win, feature_recon, param)
    rmse_vals.q =    rmse_fct(feature_win.q, feature_recon.q, param.win_length);
    rmse_vals.dq =   rmse_fct(feature_win.dq, feature_recon.dq, param.win_length);
    rmse_vals.ddq =  rmse_fct(feature_win.ddq, feature_recon.ddq, param.win_length);
    rmse_vals.dddq = rmse_fct(feature_win.dddq, feature_recon.dddq, param.win_length);
    rmse_vals.x =    rmse_fct(feature_win.x, feature_recon.x, param.win_length);
    rmse_vals.dx =   rmse_fct(feature_win.dx, feature_recon.dx, param.win_length);
    rmse_vals.ddx =  rmse_fct(feature_win.ddx, feature_recon.ddx, param.win_length);
    rmse_vals.dddx = rmse_fct(feature_win.dddx, feature_recon.dddx, param.win_length);
    rmse_vals.tau =   rmse_fct(feature_win.tau, feature_recon.tau, param.win_length);
    rmse_vals.dtau =  rmse_fct(feature_win.dtau, feature_recon.dtau, param.win_length);
    rmse_vals.ddtau = rmse_fct(feature_win.ddtau, feature_recon.ddtau, param.win_length);
%     rmse_vals.cop =   rmse_fct(feature_win.cop, feature_recon.cop, param.win_length);
%     rmse_vals.dcop =  rmse_fct(feature_win.dcop, feature_recon.dcop, param.win_length);
%     rmse_vals.ddcop = rmse_fct(feature_win.ddcop, feature_recon.ddcop, param.win_length);
%     rmse_vals.com =   rmse_fct(feature_win.com, feature_recon.com, param.win_length);
%     rmse_vals.dcom =  rmse_fct(feature_win.dcom, feature_recon.dcom, param.win_length);
%     rmse_vals.ddcom = rmse_fct(feature_win.ddcom, feature_recon.ddcom, param.win_length);
%     rmse_vals.ep =    rmse_fct(feature_win.ep, feature_recon.ep, param.win_length);
    rmse_vals.ek =    rmse_fct(feature_win.ek, feature_recon.ek, param.win_length);
    rmse_vals.geo =   rmse_fct(feature_win.geo, feature_recon.geo, param.win_length);
    rmse_vals.en =    rmse_fct(feature_win.en, feature_recon.en, param.win_length);
    
    rmse_vals.arrayName = {'rmse(q))', 'rmse(dq)', 'rmse(ddq)', 'rmse(dddq)', 'rmse(x)', 'rmse(dx)', 'rmse(ddx)', ...
        'rmse(dddx)', 'rmse(tau)', 'rmse(dtau)', 'rmse(ddtau)', 'rmse(ek)', 'rmse(geo)', 'rmse(en)'}';
    
    rmse_vals.array = [rmse_vals.q
        rmse_vals.dq
        rmse_vals.ddq
        rmse_vals.dddq
        rmse_vals.x
        rmse_vals.dx
        rmse_vals.ddx
        rmse_vals.dddx
        rmse_vals.tau
        rmse_vals.dtau
        rmse_vals.ddtau
        rmse_vals.ek
        rmse_vals.geo
        rmse_vals.en];