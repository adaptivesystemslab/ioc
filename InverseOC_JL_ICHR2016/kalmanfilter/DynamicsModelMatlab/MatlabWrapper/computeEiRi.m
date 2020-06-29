function [eis, ris] = computeEiRi(mdl,frame)
%Compute eis and ris in BASE frame, ri is the vector from joint to frame

frame = mdl.getFrameByName(frame);
dof = numel(mdl.joints);
eis = zeros(3,dof);
ris = zeros(3,dof);

T_0ee = frame.t;
for i=1:numel(mdl.joints)
    T_0i = mdl.joints(i).frame_in.t;
    ri = T_0ee(1:3,4) - T_0i(1:3,4);
    ris(:,i) = ri;
    ei = T_0i(1:3,3);
    eis(:,i) = ei;
end
end