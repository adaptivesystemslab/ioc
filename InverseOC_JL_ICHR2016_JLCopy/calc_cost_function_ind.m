function J = calc_cost_function(function_name, feature_use, param, ind)
% given the function name of interest, generate the cost function
% some cost functions were from Berret2011, or from Panchea2015

if nargin < 4
    ind = 1:size(feature_use.q, 2);
end

% len = length(ind);
len = 1;

dq = feature_use.dq(:, ind);
ddq = feature_use.ddq(:, ind);
dddq = feature_use.dddq(:, ind);
ddx = feature_use.ddx(:, ind);
dddx = feature_use.dddx(:, ind);
tau = feature_use.tau(:, ind);
dtau = feature_use.dtau(:, ind);
ddtau = feature_use.ddtau(:, ind);
geo = feature_use.geo(:, ind);
en = feature_use.en(:, ind);
ep = feature_use.ep(:, ind);
ek = feature_use.ek(:, ind);
cop = feature_use.cop(:, ind);
dcop = feature_use.dcop(:, ind);

switch function_name
    case 'ddq'        
        J = sum(sum(ddq .^ 2)) / len;

    case 'dq+ddq'        
        J = sum(sum(dq .^ 2)) / len + sum(sum(ddq .^ 2)) / len;
        
    case 'dddq'        
        J = sum(sum(dddq .^ 2)) / len;
        
    case 'ddx'
        J = sum(sum(ddx .^ 2)) / len;
        
    case 'dddx'
        J = sum(sum(dddx .^2)) / len;
        
    case 'tau'
        J = sum(sum(tau .^ 2)) / len;
        
    case 'dtau'
        J = sum(sum(dtau .^ 2)) / len;
        
    case 'ddtau' % effort        
        J = sum(sum(ddtau .^ 2)) / len;
        
    case 'kinetic_en_sqrt' % geodesic
        J = sum(sum(geo .^ 2)) / len;
        
    case 'dq_tau' % energy
        J = sum(sum(abs(en))) / len;
        
    case 'potential_en'
        J = sum(sum(ep .^ 2)) / len;
        
    case 'kinetic_en'
        J = sum(sum(ek .^ 2)) / len;
        
    case 'cop'
        J = sum(sum(cop .^ 2)) / len;
        
    case 'dcop'
        J = sum(sum(dcop .^ 2)) / len;
end