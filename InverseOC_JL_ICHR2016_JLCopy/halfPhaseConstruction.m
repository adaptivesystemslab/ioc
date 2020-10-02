function [param_docsim, ccost_array, const_x_array, const_y_array] ...
    = halfPhaseConstruction(currFilestack, param_docsim)

% set up the construnction conditions
switch param_docsim.half_phase_construction
    case 'no'
        for ind_simCon = 1:length(currFilestack.ccost_array)
            ccost_array{ind_simCon} = currFilestack.ccost_array{ind_simCon};
            const_x_array{ind_simCon} = param_docsim.const_x;
            const_y_array{ind_simCon} = param_docsim.const_y;
        end
        
    case 'yes'
        ind_sinConMarker = 0;
        for ind_simCon = 1:length(currFilestack.ccost_array)
            ind_sinConMarker = ind_sinConMarker + 1;
            ccost_array{ind_sinConMarker} = currFilestack.ccost_array{ind_simCon};
            const_x_array{ind_sinConMarker} = param_docsim.const_x([1 3]);
            const_y_array{ind_sinConMarker} = param_docsim.const_y(:, 1:2);
            
            ind_sinConMarker = ind_sinConMarker + 1;
            ccost_array{ind_sinConMarker} = currFilestack.ccost_array{ind_simCon};
            const_x_array{ind_sinConMarker} = param_docsim.const_x([1 3]);
            const_y_array{ind_sinConMarker} = param_docsim.const_y(:, 2:3);
        end
        
    case 'double' % double the middle point in the x_const
        ind_sinConMarker = 0;
        for ind_simCon = 1:length(currFilestack.ccost_array)
            ind_sinConMarker = ind_sinConMarker + 1;
            ccost_array{ind_sinConMarker} = currFilestack.ccost_array{ind_simCon};
            new_const_x = param_docsim.const_x([1 2 2 3]);
            new_const_x(2) = new_const_x(2) - 1;
            const_x_array{ind_sinConMarker} = new_const_x;
            const_y_array{ind_sinConMarker} = param_docsim.const_y(:, [1:2 2:3]);
        end
end

param_docsim.ccost_array = ccost_array; % save it for hashing purposes