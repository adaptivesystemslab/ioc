%Copyright Barshan, Ghodsi 2009
%Paper: Supervised principal component analysis: Visualization, classification and
%regression on subspaces and submanifolds.
clear
clc
close all

load sonar;
load param_dim_sonar;

[p,n] = size(X);
% isnormal = input('Is normalizetion required?(1/0)');
isnormal = 1;
if isnormal==1
    X = (X-repmat(min(X')',1,n))./(repmat(max(X')',1,n)-repmat(min(X')',1,n));
    nan_ind = find(isnan(X)==1);
    X(nan_ind) = 0;
end

ratio = 0.7;
len_tr = round(n*ratio);
indeices = randperm(n);
X = X(:,indeices);
Y = Y(indeices);

X_tr = X(:,1:len_tr);
X_ts = X(:,len_tr + 1:end);
Y_tr = Y(1:len_tr);
Y_ts = Y(len_tr + 1:end);

d = 6;
param.ktype_y = 'delta_cls';
param.kparam_y = 1;
param.ktype_x = 'rbf';
param.kparam_x = param_dim_sonar(d);
[Ztr_SPCA U] = SPCA(X_tr,Y_tr,d,param);
[Ztr_KSPCA Beta] = KSPCA(X_tr,Y_tr,d,param);

%Testing Data Kernel Computation
Ktest = repmat(0,size(X_tr,2),size(X_ts,2));
for i=1:size(X_tr,2)
    for j=1:size(X_ts,2)
        Ktest(i,j) = kernel(param.ktype_x,X_tr(:,i),X_ts(:,j),param.kparam_x,[]);
    end
end

Zts_SPCA = U'*X_ts;
Zts_KSPCA = Beta'*Ktest;


