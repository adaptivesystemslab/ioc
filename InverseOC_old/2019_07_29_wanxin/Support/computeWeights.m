function [weights, residual] = computeWeights(Hhat, numWeights)
% Obtain weights from recovery matrix using cvx
% % %     sumConstraint = 1;
% % %     [~, c] = size(Hhat);
% % % 
% % %     %cvx specification 
% % %     cvx_begin quiet 
% % %         variable x(c)    
% % %         minimize (norm(Hhat*x,2))
% % %         subject to    
% % %         sum(x(1:numWeights))==sumConstraint; %#ok<*EQEFF> 
% % %     cvx_end 
% % % 
% % %     weights=x(1:numWeights);
% % %     if size(weights,1)==1
% % %         weights=weights'; 
% % %     end
% % %     residual=cvx_optval/norm(x);
% % % end
% % % 

    sumcons=1;
    n=size(Hhat,2);
    H=Hhat'*Hhat;
    f=zeros(n,1);
    Aeq=zeros(1,n);
    Aeq(1,1:numWeights)=1;
    LB=-inf*ones(n,1);
%     LB(1:numWeights,1)=-0.05;
    LB(1:numWeights,1)=0;
    UB=inf*ones(n,1);
    UB(1:numWeights)=1;
    options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','off');
    [x,res]=quadprog(H,f,[],[],Aeq,sumcons,LB,UB,[],options);
    weights=x(1:numWeights)/sumcons;
    if size(weights,1)==1
        weights=weights';
    end
    residual=res; 