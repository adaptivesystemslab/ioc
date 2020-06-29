function error=Compute_Error(estcosts,truecost)

l=size(estcosts,1);

error=[];
for i=1:l
    est=estcosts(i,:);
    c=dot(est,truecost)/dot(est,est);
    e=norm(c*est-truecost)/norm(truecost);
    error(end+1)=e;
end
