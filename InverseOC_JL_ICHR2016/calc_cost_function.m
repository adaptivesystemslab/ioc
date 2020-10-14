function J = calc_cost_function(function_name, feature_use, param)
% given the function name of interest, generate the cost function
% some cost functions were from Berret2011, or from Panchea2015

len = size(feature_use.q, 2);

% ddq = feature_use.ddq;
% dddq = feature_use.dddq;
% x = feature_use.x;
% dx = feature_use.dx;
% ddx = feature_use.ddx;
% dddx = feature_use.dddx;
% tau = feature_use.tau;
% dtau = feature_use.dtau;
% ddtau = feature_use.ddtau;
% geo = feature_use.geo;
% en = feature_use.en;
% % ep = feature_use.ep;
% ek = feature_use.ek;
% cop = feature_use.cop;
% dcop = feature_use.dcop;

% COM = feature_use.COM;
% dCOM = feature_use.dCOM;
% ddCOM = feature_use.ddCOM;
% dddCOM = feature_use.dddCOM;
% dCOM_mag = feature_use.dCOM_mag;
% ddCOM_mag = feature_use.ddCOM_mag;
% dddCOM_mag = feature_use.dddCOM_mag;
% COM_height = feature_use.COM_height;
% dCOM_height = feature_use.dCOM_height;
% ddCOM_height = feature_use.ddCOM_height;
% dddCOM_height = feature_use.dddCOM_height;
% dToeZ = feature_use.dToeZ;
% ddToeZ = feature_use.ddToeZ;
% ToeTargX = feature_use.ToeTargX;
% dToeTargX = feature_use.dToeTargX;
% COMToeX = feature_use.COMToeX;
% dCOMToeX = feature_use.dCOMToeX;
% COMTargX = feature_use.COMTargX;
% dCOMTargX = feature_use.dCOMTargX;
% 
% ddq_leg = feature_use.ddq_leg;
% dddq_leg = feature_use.dddq_leg;
% tau_leg = feature_use.tau_leg;
% en_leg = feature_use.en_leg;
% geo_leg = feature_use.geo_leg;
% 
% ddq_arm = feature_use.ddq_arm;
% dddq_arm = feature_use.dddq_arm;
% tau_arm = feature_use.tau_arm;
% en_arm = feature_use.en_arm;
% geo_arm = feature_use.geo_arm;
% 
% ddq_tor = feature_use.ddq_tor;
% dddq_tor = feature_use.dddq_tor;
% tau_tor = feature_use.tau_tor;
% en_tor = feature_use.en_tor;
% geo_tor = feature_use.geo_tor;
% 
% ang_mom = feature_use.ang_mom;
% dang_mom = feature_use.dang_mom;
% ddang_mom = feature_use.ddang_mom;




switch function_name
%     case 'ddq'        
%         J = cost_fct_calc_sq(ddq, len);
%         J = J * param.coeff_cf.ddq;       
%         
%     case 'dddq'        
%         J = cost_fct_calc_sq(dddq, len);
%         J = J * param.coeff_cf.dddq;    
%         
%     case 'ddx'
%         J = cost_fct_calc_sq(ddx, len);
%         J = J * param.coeff_cf.ddx;    
%         
%     case 'dddx'
%         J = cost_fct_calc_sq(dddx, len);
%         J = J * param.coeff_cf.dddx;    
%         
%     case 'tau'
%         J = cost_fct_calc_sq(tau, len);
%         J = J * param.coeff_cf.tau;    
%         
%     case 'dtau'
%         J = cost_fct_calc_sq(dtau, len);
%         J = J * param.coeff_cf.dtau;    
%         
%     case 'ddtau' % effort        
%         J = cost_fct_calc_sq(ddtau, len);
%         J = J * param.coeff_cf.ddtau;    
%         
%     case 'kinetic_en_sqrt' % geodesic
%         J = cost_fct_calc_abs(geo, len);
%         J = J * param.coeff_cf.geo;    
%         
%     case 'dq_tau' % energy
%         J = cost_fct_calc_abs(en, len);
%         J = J * param.coeff_cf.en;
%         
%     case 'potential_en'
%         J = cost_fct_calc_abs(ep, len);
%         J = J * param.coeff_cf.ep;    
%         
%     case 'kinetic_en'
%         J = cost_fct_calc_abs(ek, len);
%         J = J * param.coeff_cf.ek;    
%         
%     case 'cop'
%         J = cost_fct_calc_sq(cop, len);
%         J = J * param.coeff_cf.cop;    
%         
%     case 'dcop'
%         J = cost_fct_calc_sq(dcop, len);
%         J = J * param.coeff_cf.dcop;    
%         
%     case 'kinetic_en_sqrt_sq' % geodesic
%         J = cost_fct_calc_sq(geo, len);
%         J = J * param.coeff_cf.geo;    
%         
%     case 'dq_tau_sq' % energy
%         J = cost_fct_calc_sq(en, len);
%         J = J * param.coeff_cf.en;    
%         
%     case 'potential_en_sq'
%         J = cost_fct_calc_sq(ep, len);
%         J = J * param.coeff_cf.ep;    
%         
%     case 'kinetic_en_sq'
%         J = cost_fct_calc_sq(ek, len); 
%         J = J * param.coeff_cf.ek;
%         
%     case 'ddqddx'
%         Jddq = calc_cost_function('ddq', feature_use, param);
%         Jddx = calc_cost_function('ddx', feature_use, param);
%         J = Jddq + Jddx;
%         
%     case 'ddqddx2'
%         Jddq = calc_cost_function('ddq', feature_use, param);
%         Jddx = calc_cost_function('ddx', feature_use, param);
%         J = 0.5*Jddq + 0.5*Jddx;        
   
        
%     case 'dddCOM'
%         J = cost_fct_calc_sq(dddCOM, len);
%         J = J * param.coeff_cf.dddCOM;
%     
%     % COM magnitude trajectories        
%     case 'dCOM_mag'
%         J = cost_fct_calc_sq(dCOM_mag, len);
%         J = J * param.coeff_cf.dCOM_mag;
%         
%     case 'ddCOM_mag'
%         J = cost_fct_calc_sq(ddCOM_mag, len);
%         J = J * param.coeff_cf.ddCOM_mag;
%         
%     case 'dddCOM_mag'
%         J = cost_fct_calc_sq(dddCOM_mag, len);
%         J = J * param.coeff_cf.dddCOM_mag;
%         
%     % COM height (Global Z)
%     case 'COM_height'
%         J = cost_fct_calc_sq(COM_height, len);
%         J = J * param.coeff_cf.COM_height;
%         
%     case 'dCOM_height'
%         J = cost_fct_calc_sq(dCOM_height, len);
%         J = J * param.coeff_cf.dCOM_height;
%         
%     case 'ddCOM_height'
%         J = cost_fct_calc_sq(ddCOM_height, len);
%         J = J * param.coeff_cf.ddCOM_height;
%         
%     case 'dddCOM_height'
%         J = cost_fct_calc_sq(dddCOM_height, len);
%         J = J * param.coeff_cf.dddCOM_height;
%     
%     % toe motion
%     case 'dToeZ'
%         J = cost_fct_calc_abs(dToeZ, len);
%         J = J * param.coeff_cf.dToeZ;
%         
%     case 'ddToeZ'
%         J = cost_fct_calc_abs(ddToeZ, len);
%         J = J * param.coeff_cf.ddToeZ;
%         
%     % Distances in forward jump direction (Global X)
%     case 'ToeTargX'
%         J = cost_fct_calc_abs(ToeTargX, len);
%         J = J * param.coeff_cf.ToeTargX;
%         
%     case 'dToeTargX'
%         J = cost_fct_calc_abs(dToeTargX, len);
%         J = J * param.coeff_cf.dToeTargX;
%         
%     case 'COMToeX'
%         J = cost_fct_calc_abs(COMToeX, len);
%         J = J * param.coeff_cf.COMToeX;
%     
%     case 'dCOMToeX'
%         J = cost_fct_calc_abs(dCOMToeX, len);
%         J = J * param.coeff_cf.dCOMToeX;
%         
%     case 'COMTargX'
%         J = cost_fct_calc_abs(COMTargX, len);
%         J = J * param.coeff_cf.COMTargX;
%         
%     case 'dCOMTargX'
%         J = cost_fct_calc_abs(dCOMTargX, len);
%         J = J * param.coeff_cf.dCOMTargX;
%     
%     % joint trajectories and torques separated into upper and lower body
%     case 'ddq_leg'        
%         J = cost_fct_calc_sq(ddq_leg, len);
%         J = J * param.coeff_cf.ddq_leg;       
%         
%     case 'dddq_leg'        
%         J = cost_fct_calc_sq(dddq_leg, len);
%         J = J * param.coeff_cf.dddq_leg;    
%         
%     case 'tau_leg'        
%         J = cost_fct_calc_sq(tau_leg, len);
%         J = J * param.coeff_cf.tau_leg;       
%         
%     case 'dq_tau_leg' % energy
%         J = cost_fct_calc_abs(en_leg, len);
%         J = J * param.coeff_cf.en_leg;
%     
%     case 'kinetic_en_sqrt_leg' % geodesic
%         J = cost_fct_calc_abs(geo_leg, len);
%         J = J * param.coeff_cf.geo_leg;
%     
%     
%     
%     case 'ddq_arm'        
%         J = cost_fct_calc_sq(ddq_arm, len);
%         J = J * param.coeff_cf.ddq_arm;    
%         
%     case 'dddq_arm'        
%         J = cost_fct_calc_sq(dddq_arm, len);
%         J = J * param.coeff_cf.dddq_arm;       
%         
%     case 'tau_arm'        
%         J = cost_fct_calc_sq(tau_arm, len);
%         J = J * param.coeff_cf.tau_arm;    
%         
%     case 'dq_tau_arm' % energy
%         J = cost_fct_calc_abs(en_arm, len);
%         J = J * param.coeff_cf.en_arm;
%     
%     case 'kinetic_en_sqrt_arm' % geodesic
%         J = cost_fct_calc_abs(geo_arm, len);
%         J = J * param.coeff_cf.geo_arm;  
%     
%     
%     
%     case 'ddq_tor'        
%         J = cost_fct_calc_sq(ddq_tor, len);
%         J = J * param.coeff_cf.ddq_tor;       
%         
%     case 'dddq_tor'        
%         J = cost_fct_calc_sq(dddq_tor, len);
%         J = J * param.coeff_cf.dddq_tor;    
%         
%     case 'tau_tor'        
%         J = cost_fct_calc_sq(tau_tor, len);
%         J = J * param.coeff_cf.tau_tor;    
%     
%     case 'dq_tau_tor' % energy
%         J = cost_fct_calc_abs(en_tor, len);
%         J = J * param.coeff_cf.en_tor;
%     
%     case 'kinetic_en_sqrt_tor' % geodesic
%         J = cost_fct_calc_abs(geo_tor, len);
%         J = J * param.coeff_cf.geo_tor;  
%     
%     
%     
%     case 'ang_mom'        
%         J = cost_fct_calc_abs(ang_mom, len);
%         J = J * param.coeff_cf.ang_mom;
%     
%     case 'dang_mom'        
%         J = cost_fct_calc_abs(dang_mom, len);
%         J = J * param.coeff_cf.dang_mom;
%         
%     case 'ddang_mom'        
%         J = cost_fct_calc_abs(ddang_mom, len);
%         J = J * param.coeff_cf.ddang_mom;
%             
%     case 'quantity_motion_dx'
%         dx = feature_use.dx_all;
%         J = 0;
%         for i = 1:size(dx, 1)
%             J = J + param.weights.quantity_motion_dx(i) * ...
%                 sum(abs(dx(i, :)));
%         end
%         J = J / sum(param.weights.quantity_motion_dx);
%         J = (J .^ 2) * param.coeff_cf.quantity_motion_dx;
%         
%     case 'volume_bounding_box'
%         x = feature_use.x_max;
%         vol = x(1, :) .* x(2, :) .* x(3, :);
%         J = cost_fct_calc_sq(vol, len); 
%         J = J * param.coeff_cf.volume_bounding_box;
%         
%     case 'weight_effort'
%         dx = feature_use.dx_all;         
%         eff = dx.^2;         
%         for i = 1:size(dx, 1)             
%             eff(i, :) = eff(i, :)*param.weights.weight_effort(i);         
%         end
%         effSum = sum(eff, 1); % sum across joints
%         J = max(effSum);         
%         J = J * param.coeff_cf.weight_effort;
%         
%     case 'time_effort'
%         ddx = feature_use.ddx_all;         
%         eff = sum(abs(ddx), 2) / size(ddx, 2);         
%         effMag = param.weights.time_effort(:) .* eff(:);         
%         J = (sum(effMag)^2) * param.coeff_cf.time_effort;
%         
%     case 'space_effort'
%         x = feature_use.x_all;
%         eff = [0 0 0]';
%         for i = 2:size(x, 2)
%             eff = eff + abs(x(:, i) - x(:, i-1));
%         end
%         effMag = param.weights.space_effort(:) .* eff(:);
%         J = (sum(effMag)^2) * param.coeff_cf.space_effort;
%         
%     case 'flow_effort'
%         dddx = feature_use.dddx_all;
%         eff = sum(abs(dddx), 2) / size(dddx, 2);
%         effMag = param.weights.flow_effort(:) .* eff(:);
%         J = (sum(effMag)^2) * param.coeff_cf.flow_effort;
%         
%         
%             % COM cartesian trajectories
%     case 'com'
%         J = cost_fct_calc_sq(feature_use.com, len);
%         J = J * param.coeff_cf.com;
%         
%     case 'dcom'
%         J = cost_fct_calc_sq(feature_use.dcom, len);
%         J = J * param.coeff_cf.dcom;
%         
%     case 'ddcom'
%         J = cost_fct_calc_sq(feature_use.ddcom, len);
%         J = J * param.coeff_cf.ddcom;
% 
%     case 'x_anchor_1'
%         J = cost_fct_calc_sq(feature_use.x_anchor_1, len);
%         J = J * param.coeff_cf.x_anchor_1;
%         
%     case 'x_anchor_2'
%         J = cost_fct_calc_sq(feature_use.x_anchor_2, len);
%         J = J * param.coeff_cf.x_anchor_2;
%    
%     case 'rot_anchor_1'
%         J = cost_fct_calc_sq(feature_use.rot_anchor_1, len);
%         J = J * param.coeff_cf.rot_anchor_1;
%         
%     case 'rot_anchor_2'
%         J = cost_fct_calc_sq(feature_use.rot_anchor_2, len);
%         J = J * param.coeff_cf.rot_anchor_2;
%         
%     case 'x_displace'
%         J = cost_fct_calc_sq(feature_use.x_displace, len);
%         J = J * param.coeff_cf.x_displace;
%         
%     case 'dx_displace'
%         J = cost_fct_calc_sq(feature_use.dx_displace, len);
%         J = J * param.coeff_cf.dx_displace;
%         
%     case 'ddx_displace'
%         J = cost_fct_calc_sq(feature_use.ddx_displace, len);
%         J = J * param.coeff_cf.ddx_displace;
% 
%     case 'cartCurv'
%         J = cost_fct_calc_sq(feature_use.cartCurv, len);
%         J = J * param.coeff_cf.cartCurv;
%         
%     case 'shapeDir'
%         J = cost_fct_calc_sq(feature_use.shapeDir, len);
%         J = J * param.coeff_cf.shapeDir;
%         
%     case 'x_cartCurv'
%         J = cost_fct_calc_sq(feature_use.x_cartCurv, len);
%         J = J * param.coeff_cf.x_cartCurv;
        
    case 'half_joint_task'
        J = cost_fct_calc_sq(feature_use.half_joint_task, len);
        J = J * param.coeff_cf.half_joint_task;
        
    case 'joint_length'
        J = cost_fct_calc_sq(feature_use.joint_length, len);
        J = J * param.coeff_cf.joint_length;
        
    case 'joint_limit'
        J = cost_fct_calc_sq(feature_use.joint_limit, len);
        J = J * param.coeff_cf.joint_limit;
        
    case 'manip_rot'
        J = cost_fct_calc_sq(feature_use.manip_rot, len);
        J = J * param.coeff_cf.manip_rot;
        
    case 'manip_trans'
        J = cost_fct_calc_sq(feature_use.manip_trans, len);
        J = J * param.coeff_cf.manip_trans;
        
    case 'manipulability'
        J = cost_fct_calc_sq(feature_use.manipulability, len);
        J = J * param.coeff_cf.manipulability;
        
    case 'orientation_length'
        J = cost_fct_calc_sq(feature_use.orientation_length, len);
        J = J * param.coeff_cf.orientation_length;
        
    case 'task_length'
        J = cost_fct_calc_sq(feature_use.task_length, len);
        J = J * param.coeff_cf.task_length;
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