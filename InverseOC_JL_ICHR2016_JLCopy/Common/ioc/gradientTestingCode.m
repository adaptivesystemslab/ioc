function gradientTe_symstingCode
    % running function 
    syms x1_t1 x2_t1 x3_t1 u1_t1 u2_t1 niu1_t1 niu2_t1 niu3_t1 real % time-dep variables
    syms x1_t2 x2_t2 x3_t2 u1_t2 u2_t2 niu1_t2 niu2_t2 niu3_t2 real
    syms c1_sym c2_sym phi1_sym phi2_sym Te_sym real % constant variables
    
    N = 100; % timestamps
    h = 0.1; % offset in the variable for differentiation
    
    % generaTe_sym two verions of the time-dep variables: t1 (ie var(t)) and t2
    % (ie var(t+1)) 
    e1 = x1_t1 + Te_sym*phi1_sym*u1_t1*cos(x3_t1);
    e2 = x2_t1 + Te_sym*phi1_sym*u1_t1*sin(x3_t1);
    e3 = x3_t1 + Te_sym*phi2_sym*u2_t1;
    
    % sometimes multiple independent variables are necessary, such as when
    % our model spans multiple timepoints
    indepVars{1} = [x1_t1 x2_t1 x3_t1 u1_t1 u2_t1];
    indepVars{2} = [x1_t2 x2_t2 x3_t2 u1_t2 u2_t2];
    
    % need a numerical version as well
    indepVals{1} = ones(size(indepVars{1}));
    indepVals{2} = ones(size(indepVars{2}));
    
    % set up the cost function and the constraints for the gradient calc
    depFcts{1} =    u2_t1^2; % cost function 
    depFcts{2} = c1_sym*u1_t1^2;
    depFcts{3} = c2_sym*x3_t1^2;
%     depFcts{4} = -niu1_t1*e1 + niu1_t2*x1_t2; % physics constraints
%     depFcts{5} = -niu2_t1*e2 + niu2_t2*x2_t2;
%     depFcts{6} = -niu3_t1*e3 + niu3_t2*x3_t2;
    depFcts{4} = niu1_t2*e1 - niu1_t1*x1_t2; % physics constraints (apparently correct)
    depFcts{5} = niu2_t2*e2 - niu2_t1*x2_t2;
    depFcts{6} = niu3_t2*e3 - niu3_t1*x3_t2;
    
    % define the basic cell of the gradients
    symGradSum = calcSymGradSum(depFcts, indepVars, [])'
%     numGradSum = calcNumGradSum_sym(depFcts, indepVars, indepVals, h)'


    % compare against the hand-generad Jacobian
    i = 1; 
    c1 = rand;
    c2 = rand;
    Te = rand;
    phi1 = 1;
    phi2 = 1;
    x1(i:i+1) = rand;
    x2(i:i+1) = rand;
    x3(i:i+1) = rand;
    u1(i:i+1) = rand;
    u2(i:i+1) = rand;
    niu1(i:i+1) = rand;
    niu2(i:i+1) = rand;
    niu3(i:i+1) = rand; 
    combined = [x1; x2; x3; u1; u2; niu1; niu2; niu3];
    
    handGen =  [-niu1(i) + niu1(i+1);...
                      -niu2(i) + niu2(i+1);...
                      2*c2*x3(i)-niu3(i) - Te*u1(i)*(sin(x3(i)))*niu1(i+1) + Te*u1(i)*(cos(x3(i)))*niu2(i+1) + niu3(i+1);...
                      2*c1*u1(i) + Te*cos(x3(i))*niu1(i+1) + Te*sin(x3(i))*niu2(i+1);...
                      2*u2(i) + Te*niu3(i+1)];
    
    symGen1 = subs(symGradSum, [indepVars{1} niu1_t1 niu2_t1 niu3_t1], combined(:, 1));
    symGen2 = subs(symGen1, [indepVars{2} niu1_t2 niu2_t2 niu3_t2], combined(:, 2));
    symGen3 = subs(symGen2, [c1_sym c2_sym Te_sym phi1_sym phi2_sym], [c1 c2 Te phi1 phi2]);
    
    mean(handGen - symGen3)^2
end



