function [V,d] = eigendec(A,k,how,err)

%copyright Biggs, Ghodsi 2006
error(nargchk(3,4,nargin));

if ~exist('err','var')
    err = 0;
elseif err ~= 0
    err = 1;
end

if size(A,1) ~= size(A,2)
    error('[eigendec] ERROR: input matrix is not square\n');
end
if issparse(A)
    if err
        fprintf('[eigendec] WARNING: sparse input; converting to full\n');
    end
    A = full(A);
end

[V d] = eig(A);
d = diag(d);

if ~isreal(d)
    % Complex eigenvalues!
    n = 0;
    for i = 1:length(d)
        n = n + ~isreal(d(i));
    end
    d = real(d);
    V = real(V);
    if err
        fprintf('[eigendec] WARNING: removed %d complex eigenvalues\n', n);
    end
end
n = sum(d < 0);
if n
    % Negative eigenvalues!
    if err
        fprintf('[eigendec] WARNING: %d negative eigenvalues\n', n);
    end
end

% Sort eigenvalues in descending order
[d I] = sort(d,1,'descend');
V = V(:,I);

% If 'how' is prefixed by the letter 'n', then grab all the eigenvectors with
% nonzero eigenvalues.  We use the rank of A to determine which eigenvalues
% are zero.
if strcmp(how(1), 'n')
    [sorted_abs_d I] = sort(abs(d),1,'descend');
    r = rank(A);
    f = find(sorted_abs_d > (sorted_abs_d(1) / 1e6));
    if ~isempty(f) && f(end) < r
%         fprintf('Removing zero-looking eigenvalues that made it past the rank :o\n');
        r = f(end);
    end
    I = I(1:r);
    V = V(:,I);
    d = d(I);
    how = how(2:end);
end

if strcmp(how(1),'P')   
    d = d( d > 0);
    how = how(2:end);
elseif strcmp(how(1),'N')
    d = d( d < 0);
    how = how(2:end);
end
if strcmp(how, 'LM')
    [d,I] = sort(d,1,'descend');
elseif strcmp(how, 'SM')
    [d,I] = sort(d,1,'ascend');
elseif strcmp(how, 'LA')
    [sorted_abs_d,I] = sort(abs(d),1,'descend');
    d = d(I);
elseif strcmp(how, 'SA')
    [sorted_abs_d,I] = sort(abs(d),1,'ascend');
    d = d(I);
else
    error('Invalid option "how" to eigendec');
end

V = V(:,I);
% Retrieve merely the "top" k
mk = min(k,length(d));
d = d(1:mk);
V = V(:,1:mk);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
