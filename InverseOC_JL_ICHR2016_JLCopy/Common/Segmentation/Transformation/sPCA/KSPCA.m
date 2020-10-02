function [Z Beta D] = KSPCA(X, Y, d, param)
%Copyright Barshan, Ghodsi 2009
%Paper: Supervised principal component analysis: Visualization, classification and
%regression on subspaces and submanifolds.
% Z = SPCA(X,Y,d,param)
% Input:
%       X:  explanatory variable (pxn)
%       Y:  response variables (lxn)
%       d:  number of projection dimensions
%       param: 
%             param.ktype_y : kernel type of the response variable
%             param.kparam_y : kernel parameter of the response variable
%             param.ktype_x : kernel type of the explanatory variable
%             param.kparam_x : kernel parameter of the explanatory variable


% Output:
%       Z:  dimension reduced data (dxn)
%       Beta:  U = Phi(X)xBeta where U is the orthogonal projection matrix (pxd)
if size(X,2)~=size(Y,2)
    error('X and Y must be the same length')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computing Kernel Function of Labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[l,n] = size(Y);
L = repmat(0,n,n);
%making L full rank for classification
if strcmp(param.ktype_y,'delta_cls')
     L = L+eye(n);
     param.ktype_y = 'delta';
end
for i = 1:n
    for j = 1:n
        L(i,j) = L(i,j) + kernel(param.ktype_y,Y(:,i),Y(:,j),param.kparam_y,[]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computing Kernel Function of Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = repmat(0,n,n);
for i = 1:n
    for j = 1:n
        K(i,j) = kernel(param.ktype_x,X(:,i),X(:,j),param.kparam_x,[]);
    end
end


H = eye(n)-1/n*(ones(n,n));
tmp = H*L*H*K';
[Beta D] = eigendec(tmp,d,'LM');
Z = Beta'*K;


