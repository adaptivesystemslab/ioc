function [Z U D] = SPCA(X, Y, d, param)
%Copyright Barshan, Ghodsi 2009
%Paper: Supervised principal component analysis: Visualization, classification and
%regression on subspaces and submanifolds.

% [Z U] = SPCA(X,Y,d,param)
% Input:
%       X:  explanatory variable (pxn)
%       Y:  response variables (lxn)
%       d:  dimension of effective subspaces
%       param: 
%             param.ktype_y : kernel type of the response variable
%             param.kparam_y : kernel parameter of the response variable

% Output:
%       Z:  dimension reduced data (dxn)
%       U:  orthogonal projection matrix (pxd)

if size(X,2)~=size(Y,2)
    error('X and Y must be the same length')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computing Kernel Function of Labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[l,n] = size(Y);
L = repmat(0,n,n);
classLabels = unique(Y); 
%making L full rank for classification
if strcmp(param.ktype_y,'delta_cls')
     L = L+eye(n);
     param.ktype_y = 'delta';
end
if strcmp(param.ktype_y, 'none')
    L = eye(n); % effectively just PCA
elseif strcmp(param.ktype_y, 'delta')
%     % assuming the data is sorted, the delta kernel matrix is just a block
%     % identity matrix, and can be calculated quickly instead of having to
%     % perform the inner product of each column individually
%     for ind_classLabels = 1:length(classLabels)
%         classLabelSpan = find(Y == classLabels(ind_classLabels));
%         L(classLabelSpan, classLabelSpan) = 1;
%     end

    Y_vert = repmat(Y, length(Y), 1);
    Y_hor = repmat(Y', 1, length(Y));
    L = Y_vert == Y_hor;
    L = +L; % convert the matrix to double instead of logical
    
elseif strcmp(param.ktype_y, 'linear')
%     % assuming the data is sorted, the linear kernel matrix can be
%     % calculated by piecewise multiplication of all the individual sectors
%     for ind_classLabels = 1:length(classLabels)
%         classLabelSpan{ind_classLabels} = find(Y == classLabels(ind_classLabels));
%     end
%     
%     for ind_classX = 1:length(classLabels)
%         for ind_classY = 1:length(classLabels)
%             L(classLabelSpan{ind_classX}, classLabelSpan{ind_classY}) = classLabels(ind_classX) * classLabels(ind_classY);
%         end
%     end

    L = Y'*Y;
else
    for i = 1:n
        for j = 1:n
            L(i,j) = L(i,j) + kernel(param.ktype_y,Y(:,i),Y(:,j),param.kparam_y,[]);
        end
    end
end

L = L*param.LMatrixMultiplier + eye(n)*1;

H = eye(n)-1/n*(ones(n,n));
[p,n] = size(X);
if n>p
    tmp = X*H*L*H*X';
    [U D] = eigendec(tmp,d,'LM');
else
   [u s v] = svd(L);
   phi_Y = s^.5 * v';
   tmp = phi_Y*H*X'*X*H*phi_Y';
   [V D] = eigendec(tmp,d,'LM');
   U = X*H*phi_Y'*V*inv(diag(D)^.5);
end
Z = U'*X;


