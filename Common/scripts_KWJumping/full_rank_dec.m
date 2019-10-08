function [K L] = full_rank_dec(X)
% FRD - Full rank factorization of input matrix X. 
%
% [K L] = frd(X) will write X as the product of two matrices X = KL where both K and
% L have the same rank as X
% 
% EXAMPLE 1 : 
%
% A=rand(1000);
% [K L] = frd(A);
% sum(sum(K*L - A))
% 
% ans =
% 
%    -0.0018
%
%
% EXAMPLE 2:
%
% u=[1;2;3]
% 
% u =
% 
%      1
%      2
%      3
% 
% A = u*u'
% 
% A =
% 
%      1     2     3
%      2     4     6
%      3     6     9
% 
% [K L] = frd(A);
% K*L - A
% 
% ans =
% 
%      0     0     0
%      0     0     0
%      0     0     0
% 
% rank(K)
% 
% ans =
% 
%      1

% This version: 11/3/2009

% This code works by taking X = U D V' and re writing D = D1 * D2.
% Then X = (U*D1)(D2*V').


[U S V] = svd(X) ;

[D1 D2] = splitD(S) ;

K = U*D1 ;
L = D2*V' ;

% existence of close to 0 elements may produce K and L with rank not equal to A, we want
% to round these to 0.
% Caution: if you keep too many digits in the rounding, non zero
% elements will still be picked up. If you keep too few digits in the rounding,
% accuracy may be lost.  I chose 6 decimal places as a compromise.
K = round2(K,6);
L = round2(L,6);

% For an example of the above problem, set the digits kept to 12 in the two rounding statements, then do Example 2
% from above. Rank K will equal 2, while rank A is 1.





% ---------------------------------------
% subfunctions below
% ---------------------------------------

function [D1 D2] = splitD(D)
%splitD - take a diagonal  (but not necessarily square)  matrix and "split" it into the
%product of two other diagonal (and not necessarily square) matrices, D = D1 * D2.

[n,p] = size(D);

% create placeholders for D = D1 * D2
D1 = zeros(n,p);
D2 = zeros(p,p);

a = diag(D);

sqrt_a = a.^(.5) ;

for i = 1:length(a)
    D1(i,i) = sqrt_a(i);
    D2(i,i) = sqrt_a(i);
end

function ta = round2(X,p)
% ROUND2(X,p) - round input X to p decimal places.
% If p is omitted, default is p=2.

% This method comes from Loren Shure's column "Loren on the Art of Matlab",
% September 3rd, 2009, "Rounding Results"
% http://blogs.mathworks.com/loren/2009/09/03/rounding-results/


if nargin==1
    p=2;
end

ta = round(X .* (10^p)) ./ (10^p) ;