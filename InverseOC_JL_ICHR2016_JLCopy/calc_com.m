function calc_com(q, param)

% TODO DOESNT WORK RIGHT NOW

    lengthQ = size(q, 2);
    lengthL = size(param.L);

    
    totalMass = sum(param.M(linkCount));
    
    init_mat = repmat(eye(4),[1 1 param.NbSample]);
    
    currL = [0; param.M].*param.L; % link 2 is a
    Tr = FKM_SQUAT_7DOF(init_mat, q', currL'); % end eff
    
    for ii=2:param.NbSample
        for jj = 1:8
            com_arr(:,jj) = Tr(1:3,4,jj,ii); % get the COM 
        end
        mx(:,ii) = sum(mx_arr(:, 3:5), 2)/totalMass;
    end
end