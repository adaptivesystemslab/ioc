function gradientTestingCode
    % running function 
    syms x1_t1 x2_t1 x3_t1 u1_t1 u2_t1 niu1_t1 niu2_t1 niu3_t1 real % time-dep variables
    syms x1_t2 x2_t2 x3_t2 u1_t2 u2_t2 niu1_t2 niu2_t2 niu3_t2 real
    syms c1 c2 phi1 phi2 Te real % constant variables
    
    N = 100; % timestamps
    h = 0.1; % offset in the variable for differentiation
    
    % generate two verions of the time-dep variables: t1 (ie var(t)) and t2
    % (ie var(t+1)) 
    e1 = x1_t1 + Te*phi1*u1_t1*cos(x3_t1);
    e2 = x2_t1 + Te*phi1*u1_t1*sin(x3_t1);
    e3 = x3_t1 + Te*phi2*u2_t1;
    
    % sometimes multiple independent variables are necessary, such as when
    % our model spans multiple timepoints
    indepVars{1} = [x1_t1 x2_t1 x3_t1 u1_t1 u2_t1];
    indepVars{2} = [x1_t2 x2_t2 x3_t2 u1_t2 u2_t2];
    
    % need a numerical version as well
    indepVals{1} = ones(size(indepVars{1}));
    indepVals{2} = ones(size(indepVars{2}));
    
    % set up the cost function and the constraints for the gradient calc
    depFcts{1} =    u2_t1^2; % cost function 
    depFcts{2} = c1*u1_t1^2;
    depFcts{3} = c2*x3_t1^2;
    depFcts{4} = -niu1_t1*e1 + niu1_t2*x1_t2; % physics constraints
    depFcts{5} = -niu2_t1*e2 + niu2_t2*x2_t2;
    depFcts{6} = -niu3_t1*e3 + niu3_t2*x3_t2;
    
    % define the basic cell of the gradients
    symGradSum = calcSymGradSum(depFcts, indepVars, [])'
%     numGradSum = calcNumGradSum(depFcts, indepVars, indepVals, h)'
end



