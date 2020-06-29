function J = calc_cost_function(function_name, feature_use, param)
% given the function name of interest, generate the cost function
% some cost functions were from Berret2011, or from Panchea2015

len = size(feature_use.q, 2);

ddq = feature_use.ddq;
dddq = feature_use.dddq;
ddx = feature_use.ddx;
dddx = feature_use.dddx;
tau = feature_use.tau;
dtau = feature_use.dtau;
ddtau = feature_use.ddtau;
geo = feature_use.geo;
en = feature_use.en;
ep = feature_use.ep;
ek = feature_use.ek;
cop = feature_use.cop;
dcop = feature_use.dcop;

COM = feature_use.COM;
dCOM = feature_use.dCOM;
ddCOM = feature_use.ddCOM;
dddCOM = feature_use.dddCOM;
dCOM_mag = feature_use.dCOM_mag;
ddCOM_mag = feature_use.ddCOM_mag;
dddCOM_mag = feature_use.dddCOM_mag;
COM_height = feature_use.COM_height;
dCOM_height = feature_use.dCOM_height;
ddCOM_height = feature_use.ddCOM_height;
dddCOM_height = feature_use.dddCOM_height;
dToeZ = feature_use.dToeZ;
ddToeZ = feature_use.ddToeZ;
ToeTargX = feature_use.ToeTargX;
dToeTargX = feature_use.dToeTargX;
COMToeX = feature_use.COMToeX;
dCOMToeX = feature_use.dCOMToeX;
COMTargX = feature_use.COMTargX;
dCOMTargX = feature_use.dCOMTargX;


switch function_name
    case 'ddq'        
        J = cost_fct_calc_sq(ddq, len);
        J = J * param.coeff_cf.ddq;       
        
    case 'dddq'        
        J = cost_fct_calc_sq(dddq, len);
        J = J * param.coeff_cf.dddq;    
        
    case 'ddx'
        J = cost_fct_calc_sq(ddx, len);
        J = J * param.coeff_cf.ddx;    
        
    case 'dddx'
        J = cost_fct_calc_sq(dddx, len);
        J = J * param.coeff_cf.dddx;    
        
    case 'tau'
        J = cost_fct_calc_sq(tau, len);
        J = J * param.coeff_cf.tau;    
        
    case 'dtau'
        J = cost_fct_calc_sq(dtau, len);
        J = J * param.coeff_cf.dtau;    
        
    case 'ddtau' % effort        
        J = cost_fct_calc_sq(ddtau, len);
        J = J * param.coeff_cf.ddtau;    
        
    case 'kinetic_en_sqrt' % geodesic
        J = cost_fct_calc_abs(geo, len);
        J = J * param.coeff_cf.geo;    
        
    case 'dq_tau' % energy
        J = cost_fct_calc_abs(en, len);
        J = J * param.coeff_cf.en;
        
    case 'potential_en'
        J = cost_fct_calc_abs(ep, len);
        J = J * param.coeff_cf.ep;    
        
    case 'kinetic_en'
        J = cost_fct_calc_abs(ek, len);
        J = J * param.coeff_cf.ek;    
        
    case 'cop'
        J = cost_fct_calc_sq(cop, len);
        J = J * param.coeff_cf.cop;    
        
    case 'dcop'
        J = cost_fct_calc_sq(dcop, len);
        J = J * param.coeff_cf.dcop;    
        
    case 'kinetic_en_sqrt_sq' % geodesic
        J = cost_fct_calc_sq(geo, len);
        J = J * param.coeff_cf.geo;    
        
    case 'dq_tau_sq' % energy
        J = cost_fct_calc_sq(en, len);
        J = J * param.coeff_cf.en;    
        
    case 'potential_en_sq'
        J = cost_fct_calc_sq(ep, len);
        J = J * param.coeff_cf.ep;    
        
    case 'kinetic_en_sq'
        J = cost_fct_calc_sq(ek, len); 
        J = J * param.coeff_cf.ep;    
        
    case 'ddqddx'
        Jddq = calc_cost_function('ddq', feature_use, param);
        Jddx = calc_cost_function('ddx', feature_use, param);
        J = Jddq + Jddx;
        
    case 'ddqddx2'
        Jddq = calc_cost_function('ddq', feature_use, param);
        Jddx = calc_cost_function('ddx', feature_use, param);
        J = 0.5*Jddq + 0.5*Jddx;        
    
    % COM cartesian trajectories
    case 'COM'
        J = cost_fct_calc_sq(COM, len);
        J = J * param.coeff_cf.COM;
        
    case 'dCOM'
        J = cost_fct_calc_sq(dCOM, len);
        J = J * param.coeff_cf.dCOM;
        
    case 'ddCOM'
        J = cost_fct_calc_sq(ddCOM, len);
        J = J * param.coeff_cf.ddCOM;
        
    case 'dddCOM'
        J = cost_fct_calc_sq(dddCOM, len);
        J = J * param.coeff_cf.dddCOM;
    
    % COM magnitude trajectories        
    case 'dCOM_mag'
        J = cost_fct_calc_sq(dCOM_mag, len);
        J = J * param.coeff_cf.dCOM_mag;
        
    case 'ddCOM_mag'
        J = cost_fct_calc_sq(ddCOM_mag, len);
        J = J * param.coeff_cf.ddCOM_mag;
        
    case 'dddCOM_mag'
        J = cost_fct_calc_sq(dddCOM_mag, len);
        J = J * param.coeff_cf.dddCOM_mag;
        
    % COM height (Global Z)
    case 'COM_height'
        J = cost_fct_calc_sq(COM_height, len);
        J = J * param.coeff_cf.COM_height;
        
    case 'dCOM_height'
        J = cost_fct_calc_sq(dCOM_height, len);
        J = J * param.coeff_cf.dCOM_height;
        
    case 'ddCOM_height'
        J = cost_fct_calc_sq(ddCOM_height, len);
        J = J * param.coeff_cf.ddCOM_height;
        
    case 'dddCOM_height'
        J = cost_fct_calc_sq(dddCOM_height, len);
        J = J * param.coeff_cf.dddCOM_height;
    
    % toe motion
    case 'dToeZ'
        J = cost_fct_calc_abs(dToeZ, len);
        J = J * param.coeff_cf.dToeZ;
        
    case 'ddToeZ'
        J = cost_fct_calc_abs(ddToeZ, len);
        J = J * param.coeff_cf.ddToeZ;
        
    % Distances in forward jump direction (Global X)
    case 'ToeTargX'
        J = cost_fct_calc_abs(ToeTargX, len);
        J = J * param.coeff_cf.ToeTargX;
        
    case 'dToeTargX'
        J = cost_fct_calc_abs(dToeTargX, len);
        J = J * param.coeff_cf.dToeTargX;
        
    case 'COMToeX'
        J = cost_fct_calc_abs(COMToeX, len);
        J = J * param.coeff_cf.COMToeX;
    
    case 'dCOMToeX'
        J = cost_fct_calc_abs(dCOMToeX, len);
        J = J * param.coeff_cf.dCOMToeX;
        
    case 'COMTargX'
        J = cost_fct_calc_abs(COMTargX, len);
        J = J * param.coeff_cf.COMTargX;
        
    case 'dCOMTargX'
        J = cost_fct_calc_abs(dCOMTargX, len);
        J = J * param.coeff_cf.dCOMTargX;
        
end
end

function J = cost_fct_calc_sq(feat, len)
	J = sum(sum(feat .^ 2)) / len;
% 	J = sum(sum(abs(feat))) / len;
end

function J = cost_fct_calc_abs(feat, len)
% 	J = sum(sum(feat .^ 2)) / len;
	J = sum(sum(abs(feat))) / len;
end