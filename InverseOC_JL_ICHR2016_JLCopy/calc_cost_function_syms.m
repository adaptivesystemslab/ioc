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

switch function_name
    case 'ddq'        
        J = sumSymsSq(ddq) / len;
        
    case 'dddq'        
        J = sumSymsSq(dddq) / len;
        
    case 'ddx'
        J = sumSymsSq(ddx) / len;
        
    case 'dddx'
        J = sumSymsSq(dddx) / len;
        
    case 'tau'
        J = sumSymsSq(tau) / len;
        
    case 'dtau'
        J = sumSymsSq(dtau) / len;
        
    case 'ddtau' % effort        
        J = sumSymsSq(ddtau) / len;
        
    case 'kinetic_en_sqrt' % geodesic
        J = sumSymsSq(geo) / len;
        
    case 'dq_tau' % energy
        J = sumSymsSq(abs(en)) / len;
        
    case 'potential_en'
        J = sumSymsSq(ep) / len;
        
    case 'kinetic_en'
        J = sumSymsSq(ek) / len;
        
    case 'cop'
        J = sumSymsSq(cop) / len;
        
    case 'dcop'
        J = sumSymsSq(dcop) / len;
end
end

function sumVal = sumSymsSq(symsArray)
    sumVal = sym(0);
    
    for ind1 = 1:size(symsArray, 1)
        for ind2 = 1:size(symsArray, 2)
            sumVal = sumVal + symsArray(ind1, ind2).^2;
        end
    end
end