function [R10, R01] = points2rot(P,Q,R)
% Computes the rotation matrices matrix whose column vectors are the 
% coordinates of the (unit vectors along the) axes of frame o1x1y1 defined
% by P,Q,R expressed relative to frame o0x0y0

%Now the normal vector to the plane can be calculated through the cross
%product of vectors PQ and PR 

x = (P-Q)./repmat(sqrt(sum((P-Q).^2,2)),1,3);
n = cross(x,(R-Q));
%So any point A in this plane will sattisfy the following 
% <n1,n2,n3> . (<A1,A2,A3> - <P1,P2,P3>) = 0 since the normal vector n has
% to be perpendicular to all vectors in the plane

%Now normalize vectors we will be using and find the third axis as the
%cross product on n and (Q-P)
z = n./repmat(sqrt(sum(n.^2,2)),1,3);
y = cross(z,x);
%The find rotation matrix R10 from IMU frame into world frame
R10 = zeros(3,3,size(x,1));

R10(1,1,:) = dot(x,repmat([1 0 0],size(x,1),1),2);
R10(1,2,:) = dot(y,repmat([1 0 0],size(x,1),1),2);
R10(1,3,:) = dot(z,repmat([1 0 0],size(x,1),1),2);
R10(2,1,:) = dot(x,repmat([0 1 0],size(x,1),1),2);
R10(2,2,:) = dot(y,repmat([0 1 0],size(x,1),1),2);
R10(2,3,:) = dot(z,repmat([0 1 0],size(x,1),1),2);
R10(3,1,:) = dot(x,repmat([0 0 1],size(x,1),1),2);
R10(3,2,:) = dot(y,repmat([0 0 1],size(x,1),1),2);
R10(3,3,:) = dot(z,repmat([0 0 1],size(x,1),1),2);

R01 = zeros(3,3,size(x,1));
for i=1:size(x,1)
   R01(:,:,i) =  R10(:,:,i)';
end

end