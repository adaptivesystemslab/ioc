function features_normalized = normalize_features(features_use, param)
    % normalize
%     if isfield(param, 'coeff_q')
        q = features_use.q;
        dq = features_use.dq;
        ddq = features_use.ddq;
        dddq = features_use.dddq;
        x = features_use.x;
        dx = features_use.dx;
        ddx = features_use.ddx;
        dddx = features_use.dddx;
        tau = features_use.tau;
        dtau = features_use.dtau;
        ddtau = features_use.ddtau;
        cop = features_use.cop;
        dcop = features_use.dcop;
        ddcop = features_use.ddcop;
        com = features_use.com;
        dcom = features_use.dcom;
        ddcom = features_use.ddcom;
        ep = features_use.ep;
        ek = features_use.ek;
        geo = features_use.geo;
        en = features_use.en;
        
        COM = features_use.COM;
        dCOM = features_use.dCOM;
        ddCOM = features_use.ddCOM;
        dddCOM = features_use.dddCOM;
        dCOM_mag = features_use.dCOM_mag;
        ddCOM_mag = features_use.ddCOM_mag;
        dddCOM_mag = features_use.dddCOM_mag;
        COM_height = features_use.COM_height;
        dCOM_height = features_use.dCOM_height;
        ddCOM_height = features_use.ddCOM_height;
        dddCOM_height = features_use.dddCOM_height;
        dToeZ = features_use.dToeZ;
        ddToeZ = features_use.ddToeZ;
        ToeTargX = features_use.ToeTargX;
        dToeTargX = features_use.dToeTargX;
        COMToeX = features_use.COMToeX;
        dCOMToeX = features_use.dCOMToeX;
        COMTargX = features_use.COMTargX;
        dCOMTargX = features_use.dCOMTargX;
        
        
        q   =  q   .*repmat(param.coeff_feature.q,     1, size(q, 2));
        dq  =  dq  .*repmat(param.coeff_feature.dq,    1, size(dq, 2));
        ddq =  ddq .*repmat(param.coeff_feature.ddq,   1, size(ddq, 2));
        dddq = dddq.*repmat(param.coeff_feature.dddq,  1, size(dddq, 2));

        x   =  x   .*repmat(param.coeff_feature.x,     1, size(x, 2));
        dx  =  dx  .*repmat(param.coeff_feature.dx,    1, size(dx, 2));
        ddx =  ddx .*repmat(param.coeff_feature.ddx,   1, size(ddx, 2));
        dddx = dddx.*repmat(param.coeff_feature.dddx,  1, size(dddx, 2));

        tau =   tau  .*repmat(param.coeff_feature.tau,    1, size(tau, 2));
        dtau =  dtau .*repmat(param.coeff_feature.dtau,   1, size(dtau, 2));
        ddtau = ddtau.*repmat(param.coeff_feature.ddtau,  1, size(ddtau, 2));
        
        cop =   cop  .*repmat(param.coeff_feature.cop,    1, size(cop, 2));
        dcop =  dcop .*repmat(param.coeff_feature.dcop,   1, size(dcop, 2));
        ddcop = ddcop.*repmat(param.coeff_feature.ddcop,  1, size(ddcop, 2));
        
        com =   com  .*repmat(param.coeff_feature.com,    1, size(com, 2));
        dcom =  dcom .*repmat(param.coeff_feature.dcom,   1, size(dcom, 2));
        ddcom = ddcom.*repmat(param.coeff_feature.ddcom,  1, size(ddcom, 2));

        ep  =   ep  .*repmat(param.coeff_feature.ep,    1, size(ep, 2));
        ek  =   ek  .*repmat(param.coeff_feature.ek,    1, size(ek, 2));
        geo =   geo .*repmat(param.coeff_feature.geo,   1, size(geo, 2));
        en  =   en  .*repmat(param.coeff_feature.en,    1, size(en, 2));
        
%         COM = COM  .*repmat(param.coeff_feature.COM,    1, size(COM, 2));
%         dCOM = dCOM  .*repmat(param.coeff_feature.dCOM,    1, size(dCOM, 2));
%         ddCOM = ddCOM  .*repmat(param.coeff_feature.ddCOM,    1, size(ddCOM, 2));
%         dddCOM = dddCOM  .*repmat(param.coeff_feature.dddCOM,    1, size(dddCOM, 2));
%         dCOM_mag = dCOM_mag  .*repmat(param.coeff_feature.dCOM_mag,    1, size(dCOM_mag, 2));
%         ddCOM_mag = ddCOM_mag  .*repmat(param.coeff_feature.ddCOM_mag,    1, size(ddCOM_mag, 2));
%         dddCOM_mag = dddCOM_mag  .*repmat(param.coeff_feature.dddCOM_mag,    1, size(dddCOM_mag, 2));
%         COM_height = COM_height  .*repmat(param.coeff_feature.COM_height,    1, size(COM_height, 2));
%         dCOM_height = dCOM_height  .*repmat(param.coeff_feature.dCOM_height,    1, size(dCOM_height, 2));
%         ddCOM_height = ddCOM_height  .*repmat(param.coeff_feature.ddCOM_height,    1, size(ddCOM_height, 2));
%         dddCOM_height = dddCOM_height  .*repmat(param.coeff_feature.dddCOM_height,    1, size(dddCOM_height, 2));
%         jumpToeTargX = jumpToeTargX  .*repmat(param.coeff_feature.jumpToeTargX,    1, size(jumpToeTargX, 2));
%         jumpCOMToeX = jumpCOMToeX  .*repmat(param.coeff_feature.jumpCOMToeX,    1, size(jumpCOMToeX, 2));
%         jumpCOMTargX = jumpCOMTargX  .*repmat(param.coeff_feature.jumpCOMTargX,    1, size(jumpCOMTargX, 2));
        
        
        features_normalized.q = q;
        features_normalized.dq = dq;
        features_normalized.ddq = ddq;
        features_normalized.dddq = dddq;
        features_normalized.x = x;
        features_normalized.dx = dx;
        features_normalized.ddx = ddx;
        features_normalized.dddx = dddx;
        features_normalized.tau = tau;
        features_normalized.dtau = dtau;
        features_normalized.ddtau = ddtau;
        features_normalized.cop = cop;
        features_normalized.dcop = dcop;
        features_normalized.ddcop = ddcop;
        features_normalized.com = com;
        features_normalized.dcom = dcom;
        features_normalized.ddcom = ddcom;
        features_normalized.ep = ep;
        features_normalized.ek = ek;
        features_normalized.geo = geo;
        features_normalized.en = en;
        
        features_normalized.COM = COM;
        features_normalized.dCOM = dCOM;
        features_normalized.ddCOM = ddCOM;
        features_normalized.dddCOM = dddCOM;
        features_normalized.dCOM_mag = dCOM_mag;
        features_normalized.ddCOM_mag = ddCOM_mag;
        features_normalized.dddCOM_mag = dddCOM_mag;
        features_normalized.COM_height = COM_height;
        features_normalized.dCOM_height = dCOM_height;
        features_normalized.ddCOM_height = ddCOM_height;
        features_normalized.dddCOM_height = dddCOM_height;
        features_normalized.dToeZ = dToeZ;
        features_normalized.ddToeZ = ddToeZ;
        features_normalized.ToeTargX = ToeTargX;
        features_normalized.dToeTargX = dToeTargX;
        features_normalized.COMToeX = COMToeX;
        features_normalized.dCOMToeX = dCOMToeX;
        features_normalized.COMTargX = COMTargX;
        features_normalized.dCOMTargX = dCOMTargX;
%     end
end