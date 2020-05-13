function J = analyticJacobian(angleCount, dx, dy, dz)
    % takes in an angleCount and the three position equations to compute
    % the linear analytic jacobian. Expects all angles to be in q
    % ie q1, q2, q3
    
    
    
    for i = 1:angleCount
        angleVal = ['q', num2str(i)]; % determine the current angle
        eval(['syms ', angleVal]); % create the symbolic entity
        
        J(1, i) = eval(['diff(dx, ', angleVal, ');']);
        J(2, i) = eval(['diff(dy, ', angleVal, ');']);
        J(3, i) = eval(['diff(dz, ', angleVal, ');']);
    end
end