function x = pagerank(U,G,p)
% PAGERANK  Google's PageRank
% x = pagerank(U,G,p) uses the URLs and adjacency matrix produced by SURFER,
% together with a damping factory p, (default is .85), to compute and plot
% a bar graph of page rank, and print the dominant URLs in page rank order.
% See also SURFER, SPY.

% Eliminate any self-referential links

G = G - diag(diag(G));
  
% c = out-degree, r = in-degree

[n,n] = size(G);
c = sum(G,1);
r = sum(G,2);

% Any node with no out links is temporarily linked to all nodes

for j = find(c==0);
   G(:,j) = 1;
end
   
% Solve (I - p*G*D)x = delta*e.

if nargin < 3, p = .85; end
delta = (1-p)/n;
e = ones(n,1);
D = spdiags(1./sum(G)',0,n,n);
I = speye(n,n);
x = (I - p*G*D)\(delta*e);

% Bar graph of page rank.

shg
bar(x)
title('Page Rank')

% Print URLs in page rank order.

[ignore,q] = sort(-x);
disp('     page-rank  in  out  url')
k = 1;
while (k <= n) & (x(q(k)) >= .005)
   j = q(k);
   disp(sprintf(' %3.0f %8.4f %4.0f %4.0f  %s', ...
      j,x(j),r(j),c(j),U{j}))
   k = k+1;
end
